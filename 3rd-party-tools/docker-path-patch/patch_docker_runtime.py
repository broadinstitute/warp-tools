""" 
This program will parse a WDL (using miniwdl) to find instances where the 'docker' runtime parameter is hard-coded as a string URL or specified
as a variable with a default in the caller, and replace those instances with variables that are passed up to the root workflow. This allows 
the user to easily change docker locations without having to modify the WDL. For example, changing container registries for different clouds.

Note that this code use a primitive approach to picking docker variable names. It is strongly recommended that you review the data frame output
and provide hints for better variable names. 

Note, that the errors produced about re-used docker variable names should be taken seriously as they will cause the wrong docker images to be used.

For example, 
        ERROR    Variable name 'docker' is used by multiple tasks with different docker urls: 
            ['us.gcr.io/broad-gotc-prod/dragmap:1.1.2-1.2.1-2.26.10-1.11-1643839530'
            'us.gcr.io/broad-dsp-gcr-public/base/python:3.9-debian'
            'us.gcr.io/broad-gatk/gatk:4.3.0.0']
Since these variables are using wildly differently images, it wil cause problems.

Example usage:
    python3 patch_docker_runtime.py -i "https://raw.githubusercontent.com/broadinstitute/warp/develop/pipelines/broad/dna_seq/germline/single_sample/wgs/WholeGenomeGermlineSingleSample.wdl" -o ./wgs_patch/
    python3 patch_docker_runtime.py -i "/warp/pipelines/broad/dna_seq/germline/single_sample/exome/ExomeGermlineSingleSample.wdl" -d ./hint_task_override.json -o ./exome_2_4_6_patch/
    python3 patch_docker_runtime.py -i "https://raw.githubusercontent.com/broadinstitute/gatk-sv/main/wdl/GATKSVPipelineSingleSample.wdl" -o ./gatk_sv_patch/

Hint file format:
    A single JSON struct or a list of structs that provide var_name (required) and one of digest or url. For example:
[
{ "var_name": "picard_cloud_2_23_8_docker", "digest": "None", "url": "us.gcr.io/broad-gotc-prod/picard-cloud:2.23.8" },
{ "var_name": "picard_cloud_override_docker", "digest": "mho4pzthleu5uvou2vxu3djv5eb2sgf7", "url": null },
{ "var_name": "picard_cloud_override_docker", "digest": "5balq6oq23lmha27rhncb6hx56x7yrc5", "url": null }
]

The output directory will contain:
    docker_patch.diff - A patch file that can be applied to the WDL to replace the docker runtime sections with variables
    apply_patch.sh - A bash script that will copy and apply the patch to the WDL files referenced by the root workflow
    patch_dataframe.tsv - A TSV file with the docker variable names and the docker image names, easily viewable with tidy-viewer (brew install tidy-viewer)
"""
import argparse
import json
import logging
import os
import sys
import time

import humanize
import ipdb
import pandas as pd
import WDL
from rich.pretty import pprint
from WDL.CLI import make_read_source
from WDL_Runtime_Docker_Patch import WDL_Runtime_Docker_Patch
from WDL_Task_Tree import WDL_Task_Tree, count_nodes, create_hierarchy

try:
    if get_ipython() is not None:
        import nest_asyncio
        nest_asyncio.apply()
        # Don't use the rich logger in Jupyter as this messes up Jupyter output:
        log = logging.getLogger()
        FORMAT = "%(message)s"
        logging.basicConfig(level=logging.INFO, format=FORMAT, datefmt="[%X]")
except:
    from rich.logging import RichHandler
    FORMAT = "%(message)s"
    logging.basicConfig(level="INFO", format=FORMAT, datefmt="[%X]", handlers=[RichHandler(show_path=False)])
    log = logging.getLogger("rich")


parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-i', '--input', type=str, dest='input_file', required=True,
                    help='Path to a local root WDL file to process or a http URL to a WDL file to process, note SAS tokens for imports are not currently supported')
parser.add_argument('-o', '--output', type=str, dest='output_dir', required=False, default=None,
                    help='Specify a output directory to write the "docker_patch.diff" file along with a script to perform the patch operations.')
parser.add_argument('-d', '--docker-hint', dest='docker_hint_filename', required=False, default=None,
                    help='A path to a file containing a JSON struct describing the docker image to use for each task. This is used to override default docker_naming. JSON format is "{\'var_name\': \'docker_image\', \'url\': \'docker_image_url\', \'digest\': \'wdl_task_digest\'}". Specifying var_name without a url or digest blacklists a name, for example by default "docker" is blocked as it is bad variable.')


def main():
    args = parser.parse_args()
    log.info('Calling Arguments:')  # Note: args uses the internal variable names, so we use sys instead
    log.info(" ".join("\"" + arg + "\"" if " " in arg else arg for arg in sys.argv))

    input_file, output_dir, docker_hints = setup_inputs(args)
    # Parse the root WDL (which can be either a URL or a local file) and all imports
    print("\n")
    log.info("---------------------------------")    
    log.info(f"Scanning WDL with miniWDL, this may take some time depending on workflow complexity: \"{input_file}\"")
    start_time = time.time()
    parent_wdl = WDL.load(input_file, read_source=make_read_source(False), check_quant=False)
    end_time = time.time()
    log.info(f"Successfully parsed WDL: \"{input_file}\" in {humanize.precisedelta(end_time - start_time, minimum_unit='seconds')}")

    # Create a tree structure from the WDL AST
    tree = []
    tree_parser = WDL_Task_Tree(tree)
    tree_parser.parse_root_workflow(parent_wdl.workflow)

    # Create a hierarchy from the tree
    hierarchy = create_hierarchy(tree_parser.tree)
    # log.info("Unique tree paths in WDL:")
    # pprint(tree_parser.level.get_all_paths())
    log.info(f"Counts for parsed node types: {tree_parser.stats}")
    log.info(f"Number of input nodes in the parsed tree: {len(tree_parser.tree)}")
    log.info(f"Number of nodes placed into hierarchy: {count_nodes(hierarchy[0])}")

    # Scan the WDL for docker runtime sections and create a set of patches:
    scan_docker = WDL_Runtime_Docker_Patch()
    scan_docker.scan(tree_parser, parent_wdl, docker_var_map=docker_hints)

    # Report on the patch creation:
    all_wdl_files = scan_docker.recursive_scan_imports(parent_wdl)
    wdl_files_with_patches = {k: v for k, v in scan_docker.modified_docs.items() if v.create_patch()}
    print("\n")
    log.info("---------------------------------")    
    log.info("Patch Summary:")
    log.info(f"# There are {len(all_wdl_files)} unique WDL document (including all imports)")
    log.info(f"# There are {len(wdl_files_with_patches)} modified documents")
    for name in all_wdl_files:
        if name in wdl_files_with_patches:
            log.info(f"{name} - modified")
        else:
            log.info(f"{name} - NOT modified")

    if not output_dir:
        print("\n")
        log.info("---------------------------------")
        log.info("No output directory specified, printing patch to console")
        scan_docker.print_patch()
    else:
        output_patch(output_dir, parent_wdl, scan_docker)
    log.info("Done, exiting")
    exit(0)


def output_patch(output_dir, parent_wdl, scan_docker):
    print("\n")
    log.info("---------------------------------")
    log.info(f"Writing patch script to directory: {output_dir}")

    # Write the patch to a file:
    patches, files_patch = scan_docker.get_patch()
    docker_patch_filename = 'docker_patch.diff'
    patch_filename = os.path.join(output_dir, docker_patch_filename)
    log.info(f"Writing patch to file: {patch_filename}")
    with open(patch_filename, 'w') as f:
        f.write(patches)

    # Report on patch size per files:
    common_prefix = os.path.commonprefix(list(files_patch.keys()))
    stats = {}
    for name, patch in files_patch.items():
        changed_line = 0
        num_patches = 0
        for line in patch:
            if line.startswith('+'):
                changed_line += 1
            if line.startswith('@@'):
                num_patches += 1
        clean_name = str(name).replace(common_prefix, '')
        stats[clean_name] = {'Lines Patched': changed_line, 'Number of patches': num_patches}
    df = pd.DataFrame.from_dict(stats, orient='index')
    df.index.name = "--WDL file path--"
    df.sort_values(by=['Lines Patched'], ascending=False, inplace=True)
    pd.set_option('display.max_colwidth', None)
    log.info("Patch stats per file:")
    pprint(df)

    # Write a script to apply the patch:
    script_filename = os.path.join(output_dir, 'apply_patch.sh')
    log.info(f"Writing script to apply patch to file: {script_filename}")
    log.info(f"You may run this by going to the output directory ({output_dir}) and running: ./apply_patch.sh")
    script_start = """#!/bin/bash

# This script will apply the docker patch to the WDL files in this directory
# set -e will cause the script to exit if any command fails
set -e

echo "This script will apply the docker patch to the WDL files in this directory"
echo "Original WDL: \"{source_wdl}\""


# Check if the patch has already been applied:
if [ -f ./docker_patched_files/patched ]; then
    echo "Patch has already been applied, exiting"
    exit 0
fi

# Create a directory to store the patched files:
mkdir -p ./docker_patched_files/
mkdir -p ./source_wdl_files/

# Copy the original files to the patched directory:
echo "Copying original files to ./docker_patched_files/"
echo "Unmodified/unpatched files will also be copied to ./source_wdl_files/"
echo "Note by running diff on these two directories you can see the changes made by the patch"
echo "Such as: diff -ur ./docker_patched_files/ ./source_wdl_files/"
echo ""
echo "Copying files..."
""".format(source_wdl=parent_wdl.pos.abspath)

    mkdir_ops = ""
    download_ops = ""
    copy_ops = ""
    for name in scan_docker.wdl_file_list:
        if name.startswith("http:") or name.startswith("https:"):
            # We have to use curl to download the file:
            source = name
            name = name.replace("http:/", '')
            name = name.replace("https:/", '')
            
            if not name.startswith('/'):
                name = '/' + name
            dest = name.replace(scan_docker.abspath_root, './docker_patched_files')
            clean_dest = dest.replace('./docker_patched_files', '')
            mkdir_ops += f'mkdir -p "$(dirname "{dest}")"\n'
            download_ops += f'curl -s "{source}" > "{dest}"\n'
            mkdir_ops += f'mkdir -p "$(dirname "./source_wdl_files{clean_dest}")"\n'
            copy_ops += f'cp "{dest}" "./source_wdl_files{clean_dest}"\n'
        else:
            source = name
            if not name.startswith('/'):
                name = '/' + name
            dest = name.replace(scan_docker.abspath_root, './docker_patched_files')
            clean_dest = dest.replace('./docker_patched_files', '')
            mkdir_ops += f'mkdir -p "$(dirname "{dest}")"\n'
            copy_ops += f'cp "{source}" "{dest}"\n'
            mkdir_ops += f'mkdir -p "$(dirname "./source_wdl_files{clean_dest}")"\n'
            copy_ops += f'cp "{dest}" "./source_wdl_files{clean_dest}"\n'
    
    source_wdl_filename = parent_wdl.pos.abspath
    if source_wdl_filename.startswith("http:") or source_wdl_filename.startswith("https:"):
        source_wdl_filename = source_wdl_filename.replace("http:/", '')
        source_wdl_filename = source_wdl_filename.replace("https:/", '')
        source_wdl_filename = source_wdl_filename.replace(scan_docker.abspath_root, './docker_patched_files')
    else:
        source_wdl_filename = source_wdl_filename.replace(scan_docker.abspath_root, './docker_patched_files')

    script_end = """
echo "Done copying files"
echo ""

# Apply the patch:
echo ""
echo "Applying patch"
patch -d ./docker_patched_files/ -p1 < {patch_filename}
echo "Done applying patch"

# Create a file to indicate the patch has been applied:
touch ./docker_patched_files/patched

echo ""
echo "You may validate the patched WDL by using miniwdl or womtool"
echo "For example: miniwdl check --suppress CommandShellCheck --no-quant-check {source_wdl_filename}"
echo ""

echo ""
echo "Checking with miniwdl"
if miniwdl check --suppress CommandShellCheck --no-quant-check  "{source_wdl_filename}" >/dev/null
then
    echo "Success!"
else
    echo "Error!!"
fi

echo ""
echo "Done, exiting"
exit 0
""".format(patch_filename=docker_patch_filename,source_wdl_filename=source_wdl_filename)
    with open(script_filename, 'w') as f:
        f.write(script_start)
        f.write(mkdir_ops)
        f.write(download_ops)
        f.write(copy_ops)
        f.write(script_end)

    # Write a TSV file with the docker variable names and the docker image names:
    # easily viewable with tidy-viewer (brew install tidy-viewer)
    # tv patch_dataframe.tsv -e -u 90
    dataframe_filename = os.path.join(output_dir, 'patch_dataframe.tsv')
    scan_docker.var_name_df.to_csv(dataframe_filename, sep='\t', index=False)

    # ipdb.set_trace()


def setup_inputs(args):
    input_file = args.input_file
    output_dir = args.output_dir

    # Make sure the output directory exists:
    if output_dir and os.path.exists(output_dir) is False:
        log.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    # Load the docker_hints JSON file if it was supplied and exists:
    docker_hint_filename = args.docker_hint_filename
    if docker_hint_filename and os.path.exists(docker_hint_filename):
        print("\n")
        log.info("---------------------------------")
        log.info(f'Loading docker hints from file: "{docker_hint_filename}"')
        with open(docker_hint_filename) as f:
            docker_hints = json.load(f)
        if isinstance(docker_hints, dict):
            docker_hints = [docker_hints]
    elif docker_hint_filename and os.path.exists(docker_hint_filename) is False:
        log.error(f"Specified docker hint file does not exist: {docker_hint_filename}")
        sys.exit(1)
    else:
        docker_hints = None

    # Perform some simple validation on the JSON supplied:
    if docker_hints:
        # Convert any None strings to actual None:
        docker_hints = [{k: None if v == 'None' else v for k, v in json_obj.items()} for json_obj in docker_hints]

        # For a given row if digest and url are both specified, warn that only the digest will be used:
        docker_df = pd.DataFrame(docker_hints)
        docker_df['digest_and_url'] = docker_df.apply(lambda row: True if row['digest'] and row['url'] else False, axis=1)
        if docker_df['digest_and_url'].any():
            log.warning("WARNING: For the following rows, the digest will be used and the url will be ignored:")
            pprint(docker_df[docker_df['digest_and_url'] == True])
            print("\n")
        docker_df.drop('digest_and_url', axis=1, inplace=True)
        log.info("Docker variable name overrides loaded from file:")
        pprint(docker_df)

    return input_file, output_dir, docker_hints


if __name__ == '__main__':
    main()
