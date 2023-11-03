import difflib
import datetime
import requests
import os
import logging
import WDL
import pandas as pd

from Level import Level
log = logging.getLogger()


# TODO: WDL does not actually require that all inputs and parameters be defined as the only thing on a newline (whitespace doesn't matter), 
# both the patch class and the docker runtime variable patcher this code assumes that everything is on a newline. To fix this, patches must 
# use the line and column/end_column ranges to constrain changes.
class WDL_Document_Patch():
    def __init__(self, wdl_data, abspath_root=None):
        self.original_lines = [x + '\n' for x in wdl_data.source_lines]
        self.modified_lines = self.original_lines.copy()
        self.lines_to_drop = []   # List of lines to empty before diffing (to avoid breaking line numbers)
        # Add a blank line to indicate the end of the file (otherwise difflib produces extra changes at the end of a file)
        self.doc = wdl_data
        self.original_lines.append('')
        self.modified_lines.append('')
        # Store the path for the file:
        if not abspath_root:
            self.original_path = wdl_data.pos.uri
            self.original_abspath = wdl_data.pos.abspath
            self.modified_path = wdl_data.pos.uri
        else:
            self.original_path = wdl_data.pos.abspath.replace(abspath_root, '.')
            self.original_abspath = wdl_data.pos.abspath
            self.modified_path = wdl_data.pos.abspath.replace(abspath_root, '.')
        return

    def drop_line_pairs(self, modified) -> list:
        # Drop lines in reverse order to avoid breaking our indices during processing:
        drop_list = sorted(self.lines_to_drop, key=lambda x: x[1], reverse=True)
        for line in drop_list:
            modified[line[0]:line[1]] = ''
        return modified

    def expand_new_lines(self, modified) -> list:
        # If any element in the list of lines contains more than one new line, split it into enough lines
        # such that there's only a single new line per element. Done in reverse order to avoid indexing issues
        expand_lines = [idx for idx, line in enumerate(modified) if line.count('\n') > 1]
        expand_lines = sorted(expand_lines, reverse=True)

        for line in expand_lines:
            # Split by new line
            new_lines = modified[line].split('\n')
            # if len(new_lines) > 2:
            #     # Remove duplicate new lines
            #     new_lines = list(set(new_lines))
            # Remove the trailing new line
            new_lines.remove('')
            # Take the first part of the file and the last part (minus the line we expanded) and insert a new set of lines in the middle
            modified = modified[0:line] + list(map(lambda x: x + '\n', new_lines)) + modified[line + 1:]
        return modified
    
    def create_patch(self, context_lines=3):
        modified = self.modified_lines.copy()
        if self.lines_to_drop:  # If we have lines to drop, we need to drop them before diffing
            modified = self.drop_line_pairs(modified)
        # We need to expand any lines which contain more than one new line before diffing
        modified = self.expand_new_lines(modified)

        # Attempt to get the modified time of the original file
        if self.original_abspath.startswith('http://') or self.original_abspath.startswith('https://'):
            try:  # See if we can get the last modified time from the server
                original_mtime = requests.head(self.original_abspath).headers['Last-Modified']
            except:
                original_mtime = datetime.datetime.now().astimezone().isoformat()
        else:
            original_mtime = datetime.datetime.fromtimestamp(os.path.getmtime(self.original_abspath), datetime.timezone.utc).astimezone().isoformat()
        modified_mtime = datetime.datetime.now().astimezone().isoformat()

        # Now we can diff the two lists of lines
        diff = difflib.unified_diff(self.original_lines, modified,
                                    fromfile=self.original_path, tofile=self.modified_path,
                                    fromfiledate=original_mtime, tofiledate=modified_mtime, n=context_lines)
        diff = [x for x in diff]
        return diff

def skip_processed_files(func):
    """Decorator to skip documents which have already been encountered as imports can occur many times
    2nd parameter must be a WDL object with pos.abspath (typically this must be a WDL.Tree.Document)
    """
    processed_files = set()
    def wrapper(*args, **kwargs):
        file_path = args[2].pos.abspath
        if file_path in processed_files:
            return None
        result = func(*args, **kwargs)
        processed_files.add(file_path)
        return result
    return wrapper

def show_dataframe(df):
    try:
        from IPython.display import display
        display(df)
    except:
        print(df)

class WDL_Runtime_Docker_Patch():
    def __init__(self):
        self.modified_docs = {}  # Index here is the abspath
        self.modified_tasks = {}  # Index here is {document.pos.abspath}_{task.name}_{task.pos.line}
        self.modified_callers = {}  # Index here is the caller's digest
        self.abspath_root = None  # This is the location of the repo (the most common abspath root if not provided)
        self.found_docker_images = {}  # Index is task digest, results from initial scan
        self.var_name_df = pd.DataFrame()  # Final naming dataframe
        self.used_docker_var_names = {}     # variable name is the key with the value being the docker url if available, otherwise the specific task digest
                                            # or None (indicating that no one else can use this variable name)
        self.wdl_file_list = []  # List of all WDL files that are used by this WDL (not all files will need to be patched)

    def scan(self, tree_parser, wdl_doc, docker_var_map=[{'var_name':'docker', 'url':None, 'digest':None}], abspath_root=None):
        """Main scanning function to read through a WDL document using a calling hierarchy"""
        # One time setup:
        self.init_setup(tree_parser, wdl_doc, abspath_root)

        # Populate a dataframe that tells us what docker_variable name to use when one is not provided
        print("\n")
        log.info("---------------------------------")
        log.info("Creating docker variable name map, please carefully check warnings and the final table for correctness:")
        self.init_create_docker_var_map(tree_parser, docker_var_map)

        # Now we can actually scan to create a patch:
        print("\n")
        log.info("---------------------------------")
        log.info("Performing actual patching of WDL documents:")
        self.scan_patch(tree_parser, wdl_doc, abspath_root)

        # Print the breakdown of docker variable names used:
        print("\n")
        log.info("---------------------------------")
        log.info("Final docker variable name usage and counts:")
        self.print_docker_var_name_df()
        return

    def print_patch(self):
        patch_output = self.get_patch()
        if patch_output:
            print(patch_output)

    def get_patch(self):
        patch_output = ""
        file_patch = {}
        for name, doc in self.modified_docs.items():
            patch = doc.create_patch()
            if patch:  # Avoid blank lines for empty patches
                patch_output += "".join(patch)
                file_patch[name] = patch
        return patch_output, file_patch

    def print_docker_var_name_df(self):
        docker_urls_by_var = self.var_name_df[['accepted_var_name', 'docker_url']].value_counts()
        docker_urls_by_var.index = docker_urls_by_var.index.set_names(['Docker Runtime Var Name', 'Docker URL'])
        show_dataframe(docker_urls_by_var)
        print("\n")
        log.info("Inputs that are now required to run this Workflow (you will now need to provide these as RootWorkflowName.variable inputs):")
        var_groups = self.var_name_df.groupby('accepted_var_name')
        # pd.set_option('display.max_colwidth', None)
        # pd.set_option('display.width', 1000)
        for name, group in var_groups:
            if group['docker_url'].isnull().all() and not group['var_name_override'].any():
                continue  # Skip if no docker url and no override

            unique_urls = set(group['docker_url'])
            if len(unique_urls) <= 1 or (len(unique_urls) == 2 and None in unique_urls):
                if not unique_urls:
                    print(f'    "{name}" : "DEFINE_THIS_URL_{name.toUpperCase()}",')
                else:
                    print(f'    "{name}" : "{group["docker_url"].iloc[0]}",')
            else:
                log.warning(f"Multiple docker urls found for {name}: {unique_urls},")
                print(f'    "{name}" : "{unique_urls}"')
        return

    def init_setup(self, tree_parser, wdl_doc, abspath_root):
        if abspath_root:
            # User supplied abspath_root
            self.abspath_root = abspath_root
        else:
            # Figure out the abspath root for every WDL document, by taking the common root + up 1 level:
            abspath_list = set([x['position'].abspath for x in tree_parser.tree])
            # If we have any https:// we can't do os.path.abspath, so just use the commonpath:
            if any([x.startswith('https://') for x in abspath_list]) or any([x.startswith('http://') for x in abspath_list]):
                # Remove any http:/ or https:/ from the beginning of the abspath:
                abspath_list = [x.replace('http:/', '') for x in abspath_list]
                abspath_list = [x.replace('https:/', '') for x in abspath_list]
                common_path = os.path.commonpath(abspath_list)
                self.abspath_root = common_path
            else:
                common_path = os.path.commonpath(abspath_list)
                self.abspath_root = os.path.abspath(os.path.join(common_path, '..'))  # Ascend up a level to preserve more context
        log.info(f"Common abspath root found to be: {self.abspath_root}")

        # Scan the WDL imports to build a list of files that are used by this WDL:
        self.wdl_file_list = self.recursive_scan_imports(wdl_doc)

        # Scan the WDL document for existing docker strings and initialize the wdl patching object
        self.init_scan_docs(tree_parser, wdl_doc)


    @skip_processed_files
    def init_scan_docs(self, tree_parser, wdl_doc):
        """Iterate through all the tasks in the document and find any docker String inputs and their corresponding value/variable name """
        # Ensure the document has been initialized with the correct abspath_root
        if not self.abspath_root:
            assert self.abspath_root, f"abspath_root not set for {wdl_doc.pos.abspath}"
        meta = WDL_Document_Patch(wdl_doc, self.abspath_root)
        self.modified_docs[wdl_doc.pos.abspath] = meta

        # Process the current document
        self.init_scan_doc_for_docker_images(tree_parser, wdl_doc)

        # Process all imports:
        for document in wdl_doc.imports:
            # Get the WDL.Tree.Document from the DocImport via document.doc and recursively scan from here
            self.init_scan_docs(tree_parser, document.doc)
        return

    def init_scan_doc_for_docker_images(self, tree_parser, document):
        for task in document.tasks:
            if 'docker' not in task.runtime:
                continue
            match(type(task.runtime['docker'])):
                case WDL.Expr.String:
                    # [3] - docker image is defined as a string default in the runtime attributes:
                    default_docker_image = str(task.runtime['docker'])
                    self.found_docker_images[task.digest] = self.get_dict_task_metadata(None, default_docker_image, task)
                case WDL.Expr.Get:
                    # Docker image is defined as a variable in the runtime attributes:
                    #  Example: docker = gatk_docker
                    # Examine the variable being referenced:
                    var = task.runtime['docker'].expr.referee
                    if type(var) == WDL.Tree.Decl:
                        if not var.expr:
                            # [1] - Docker image is defined as an input without a default (which is what we want)
                            self.found_docker_images[task.digest] = self.get_dict_task_metadata(task.runtime['docker'].expr.referee.name, None, task)
                            continue
                        elif type(var.expr) == WDL.Expr.String:
                            # [2] - Docker image is defined as an task input with a default in the task:
                            docker_var_name = task.runtime['docker'].expr.referee.name
                            docker_input, docker_variable_name, default_docker_image = self.get_task_input(task, docker_var_name)
                            self.found_docker_images[task.digest] = self.get_dict_task_metadata(docker_variable_name, default_docker_image, task)
                        else:
                            assert False, f"Unknown type of expression: {type(var.expr)}, {task}"
                    else:
                        assert False, f"Variable is not a decl: {str(var)}, {task}"
                case _:
                    assert False, f"Unknown type of docker image variable: {type(task.runtime['docker'])}, {task}"

    def extract_docker_image_name(self,url):
        if url:
            return url.split('/')[-1].split(':')[0].split('@')[0]
        else:
            return None

    def get_task_input(self, task, docker_var_name):
        # Get the corresponding docker task input:
        docker_input = [x for x in task.inputs if x.name == docker_var_name]
        assert len(docker_input) == 1, f"Expected to find a single input with name {docker_var_name}, found {len(docker_input)}"
        docker_input = docker_input[0]

        # The value of the docker string should be a child of the task input:
        children = [x for x in docker_input.children]
        assert len(children) == 1, f"Expected to find a single child of the docker input, found {len(children)}"
        default_docker_image = str(children[0])
        docker_variable_name = str(docker_input.name)
        return docker_input, docker_variable_name, default_docker_image

    def get_dict_task_metadata(self, var_name, docker_url, task):
        return {'var_name': var_name, 'docker_url': docker_url, 'abspath': task.pos.abspath, 'lines': f"{task.pos.line}-{task.pos.end_line}", "task_name": task.name}

    def init_create_docker_var_map(self, tree_parser, docker_var_map):
        # Use the results from the initial scan + user input to come up with docker variable names:
        df = pd.DataFrame.from_dict(self.found_docker_images, orient='index').rename_axis('task_digest')
        # Remove self.abspath_root from the abspath column:
        # df['abspath'] = df['abspath'].str.replace(self.abspath_root, '.', regex=False)
        # Clean up quoting around docker_urls:
        df['docker_url'] = df['docker_url'].apply(lambda x: x.strip('"') if type(x) is str else x)
        patch_df = df[df['docker_url'].isnull()]  # List of docker urls some with variable names
        input_df = df[df['docker_url'].notnull()]  # list of docker variables without urls

        # Cases we care about:
        # 1. User supplied docker variable name + task digest (docker url is ignored): use this name only for specific tasks which match this digest
        # 2. User supplied docker variable name + docker url: use this name for anything which matches this url (digest must be None)
        # 3. User supplied docker variable names: don't use this name for anything else
        # 4. Scanned docker names without a docker url: blacklist this variable name to avoid confusion
        # 5. Scanned docker names with a docker url: use this name for anything which matches this url (unless it conflicts with user supplied names or is just 'docker')
        # Now using the above cases we can come with names for the unassigned docker variables, we'll use the docker image name or if this conflicts with another image add _1, _2, or the task digest if we can't get the docker url
        # naming_df  = df[df['var_name'].isnull()].copy()
        naming_df = df.copy()
        naming_df['accepted_var_name'] = None
        naming_df['var_name_source'] = None

        # Handle user supplied docker parameters (these rules are applied without checking for conflicts with other docker var names)
        for mapping in docker_var_map:
            if mapping['digest'] and mapping['var_name']:
                self.used_docker_var_names[mapping['var_name']] = mapping['digest']  # Case 1, Name is only used for this specific task
                naming_df.loc[naming_df.index.isin([mapping['digest']]), ['accepted_var_name', 'var_name_source']] = [mapping['var_name'], "user supplied by digest"]
            elif mapping['url'] and mapping['var_name'] and not mapping['digest']:
                self.used_docker_var_names[mapping['var_name']] = mapping['url']  # Case 2, Name may be used for the same docker url
                naming_df.loc[naming_df['docker_url'].isin([mapping['url']]), ['accepted_var_name', 'var_name_source']] = [mapping['var_name'], "user supplied by docker url"]
            elif mapping['var_name']:
                self.used_docker_var_names[mapping['var_name']] = None  # Case 3, Name is not allowed to be used
            else:
                log.warning(f"Unmatched docker variable mapping, ignoring: {mapping}")

        # Warn if we found inconsistent use of a variable name (that is, same variable name but different docker images):
        var_name_groups = input_df.groupby('var_name').agg({'docker_url': 'nunique'})
        for var_name, count in var_name_groups['docker_url'].items():
            if count > 1:
                log.error(f"Variable name '{var_name}' is used by multiple tasks with different docker urls: {input_df[input_df['var_name'] == var_name]['docker_url'].unique()}")

        # Now handle case 4 existing docker variable names:
        if not patch_df.empty:
            var_groups = patch_df.groupby('var_name')
            for used_var_name, group in var_groups:
                if used_var_name not in self.used_docker_var_names:
                    self.used_docker_var_names[used_var_name] = None  # Case 4, Name is not allowed to be used

        # Handle scanned docker images for cases 5:
        if not input_df.empty and input_df['docker_url'].notnull().any():
            url_groups = input_df.groupby('docker_url')
            for docker_url, group in url_groups:
                # For every docker url check to see if we have been given a variable name for this url:
                accepted_var_name = None
                for key,val in self.used_docker_var_names.items():
                    if val == docker_url:
                        accepted_var_name = key
                        var_source = "user supplied by docker url"

                # For every value where we have the same docker_url and at least one variable, attempt to use that variable name (if we are allowed to)
                if not accepted_var_name and group['var_name'].notnull().any():
                    # Accept a variable name that does not conflict with user supplied parameters:
                    var_name_counts = group['var_name'].value_counts()
                    accepted_var_name = None
                    for docker_var_name, count in var_name_counts.items():
                        if docker_var_name not in self.used_docker_var_names and docker_var_name != 'docker':
                            accepted_var_name = docker_var_name
                            var_source = "scanned, docker url, common variable name"
                            break

                # No valid variable name found, try using the image name:
                if not accepted_var_name:
                    docker_image_name = self.extract_docker_image_name(docker_url)
                    if docker_image_name:
                        accepted_var_name = self.get_docker_var_name(docker_image_name)
                        var_source = "scanned, docker url, docker image name"

                #  If we find a variable name that doesn't conflict with user supplied parameters, use it for this url:
                if accepted_var_name and accepted_var_name not in self.used_docker_var_names:
                    naming_df.loc[naming_df['docker_url'] == docker_url, ['accepted_var_name', 'var_name_source']] = [accepted_var_name, var_source]
                    # self.docker_url_to_docker_var_name[docker_url] = accepted_var_name
                    self.used_docker_var_names[accepted_var_name] = docker_url

        # Final clean-up step. For all dockers without a variable name, come up with a name (this should be only be a few or to resolve conflicts)
        for idx, row in naming_df.iterrows():
            if row['accepted_var_name']:
                continue
            # If we get here we need to come up with a variable name:
            # First try to use the image name:
            docker_image_name = self.extract_docker_image_name(row['docker_url'])
            if docker_image_name:
                accepted_var_name = self.get_docker_var_name(docker_image_name)
                var_source = "scanned, docker image name"
            else:
                # If we can't find an image name, just use a generic name:
                accepted_var_name = self.get_docker_var_name(f"{idx}_docker")
                var_source = "scanned, generic"
            naming_df.loc[naming_df.index == idx, ['accepted_var_name', 'var_name_source']] = [accepted_var_name, var_source]
            self.used_docker_var_names[accepted_var_name] = None

        # Provide feedback and create a final naming dataframe for deciding on variable names:
        fixed_var_name_df = df[df['var_name'].notnull()].copy()
        provided_var_name_rows = fixed_var_name_df['var_name'].notnull()
        fixed_var_name_df.loc[provided_var_name_rows, 'accepted_var_name'] = fixed_var_name_df.loc[provided_var_name_rows, 'var_name']
        fixed_var_name_df.loc[provided_var_name_rows, 'var_name_source'] = "existing docker variable name"
        log.info(f"Found {len(fixed_var_name_df)} docker images with existing variable names:")
        show_dataframe(fixed_var_name_df)

        log.info(f"Found {len(naming_df)} docker images with & without variable names:")
        log.info("accepted_var_name will be used when no appropriate var_name is found, or if a variable must be passed beyond it's current passing scope")
        show_dataframe(naming_df)
        log.info("---------------------------------")

        # Save the final naming dataframe:
        self.var_name_df = naming_df
        self.var_name_df['var_name_override'] = False

        # For tasks where we had a variable name + default docker image, the caller might not have been providing the docker input, we need to
        # add a docker input in those case because we are removing the default docker image
        log.info(f"Checking if docker variables are correctly passed for tasks where the docker image is defined as a variable")
        for idx, row in naming_df.iterrows():
            if not row['var_name']:
                continue
            log.info(f"   Checking {row['task_name']}: '{row['var_name']}'")
            callers = [x for x in tree_parser.tree_map.keys() if x.find(idx) != -1]
            for caller_idx in callers:
                caller = tree_parser.tree_map[caller_idx]
                override = False
                if row['var_name'] not in caller.inputs:  # The caller is not providing an input so we need to pass this in
                    log.error(f"       Error, caller '{caller.name}' not providing input ({caller_idx})")
                    override = True
                else:  # We do not check the full call tree, just the immediate caller
                    log.info(f"       Caller '{caller.name}' ({caller_idx}) providing input, '{row['var_name']}'")
                if override:
                    self.var_name_df.loc[self.var_name_df.index == idx, 'var_name_override'] = True

    def get_docker_var_name(self, docker_image_name):
        accepted_var_name = docker_image_name.replace('-', '_') + '_docker'
        existing_var_count = 1
        while accepted_var_name in self.used_docker_var_names:
            accepted_var_name = docker_image_name.replace('-', '_') + f'_docker_{existing_var_count}'
            existing_var_count += 1
        return accepted_var_name

    def get_common_abspath_root(self, wdl_doc):
        abspath_list = [wdl_doc.pos.abspath]
        abspath_list.extend(self.recursive_scan_imports(wdl_doc, True))
        common_path = os.path.commonpath(abspath_list)
        # if args.root_up_one:  # Preserve a bit more context
        common_path = os.path.abspath(os.path.join(common_path, '..'))
        self.abspath_root = common_path
        log.warning(f"Common abspath root found to be: {self.abspath_root}")

    def recursive_scan_imports(self, wdl_doc, get_abs_path=True):
        abspath_list = []
        def add_doc(wdl_doc):
            if get_abs_path:
                abspath_list.append(wdl_doc.pos.abspath)
            else:
                abspath_list.append(wdl_doc.pos.uri)

        add_doc(wdl_doc)
        for document in wdl_doc.imports:
            add_doc(document.doc)
            if document.doc.imports:
                abspath_list.extend(self.recursive_scan_imports(document.doc))
        return sorted(list(set(abspath_list)))

    def scan_patch(self, tree_parser, wdl_doc, abspath_root=None):
        # Process the current document
        self.scan_single_document(tree_parser, wdl_doc)

        # Process all imports:
        for document in wdl_doc.imports:
            # Get the WDL.Tree.Document from the DocImport via document.doc and recursively scan from here
            self.scan_patch(tree_parser, document.doc, self.abspath_root)
        return

    # TODO: remove abspath_root from here
    @skip_processed_files
    def scan_single_document(self, tree_parser, document):
        """ Given a call hierarchy of a WDL document, scan all the tasks in all imported documents and examine each docker image
            - [1] Docker images defined as inputs without a String default are left alone
                    runtime {
                        docker: gatk_docker
                        ...
            - [2] Docker images defined as inputs with a String default have their default value removed and (optionally) the call hierarchy updated
                    input {
                        String gatk_docker = "broadinstitute/gatk:
                        ...
                    runtime {
                        docker: gatk_docker
                        ...
            - [3] Docker images defined in the runtime section are changed to be an input, the input added to the task inputs, and the call hierarchy updated
                    runtime {
                        docker: "broadinstitute/gatk:

            Outputs:
                - List of patches to apply to the WDL document
                - List of modified tasks and their associated docker image defaults
            Following this scan a similar scan should be done for all callers of the tasks found to ensure defaults are not being set in callers
        """
        # log.info(f"Scanning document: {os.path.basename(document.pos.abspath)}")
        for task in document.tasks:
            if 'docker' not in task.runtime:
                continue
            match(type(task.runtime['docker'])):
                case WDL.Expr.String:
                    # [3] - docker image is defined as a string default in the runtime attributes:
                    # default_docker_image = str(task.runtime['docker'])
                    old_docker_var_name, new_docker_var_name = self.get_new_docker_var_name(task)  # Get the name for the final docker variable
                    self.patch_task_runtime_parameter(task, 'docker', new_docker_var_name)   # Fix the docker runtime section
                    self.patch_task_input(task, old_docker_var_name, new_docker_var_name)  # Fix the task inputs
                    # Now that this particular task has been patched, all callers (up the entire call tree) have to be updated:
                    self.patch_update_caller_hierarchy(tree_parser, task, old_docker_var_name, new_docker_var_name)

                    # docker_var_name, default_docker_image = self.convert_runtime_to_input(document, task)
                    # self.add_task_input(document, task, docker_var_name)
                    # self.add_input_to_hierarchy(tree_parser, task, docker_var_name, default_docker_image)
                    continue
                case WDL.Expr.Get:
                    # Docker image is defined as a variable in the runtime attributes:
                    #  Example: docker = gatk_docker
                    # Examine the variable being referenced:
                    var = task.runtime['docker'].expr.referee
                    if type(var) == WDL.Tree.Decl:
                        if not var.expr:
                            # [1] - Docker image is defined as an input without a default (which is what we want)
                            continue
                        elif type(var.expr) == WDL.Expr.String:
                            # [2] - Docker image is defined as an task input with a default in the task:
                            old_docker_var_name, new_docker_var_name = self.get_new_docker_var_name(task)  # Get the name for the final docker variable
                            self.patch_task_runtime_parameter(task, 'docker', new_docker_var_name)   # Fix the docker runtime section
                            self.patch_task_input(task, old_docker_var_name, new_docker_var_name)  # Fix the task inputs
                            # Now that this particular task has been patched, all callers (up the entire call tree) have to be updated:
                            self.patch_update_caller_hierarchy(tree_parser, task, old_docker_var_name, new_docker_var_name)
                            # if self.var_name_df.loc[task.digest, 'var_name_override']:
                            #     docker_var_name = self.var_name_df.loc[task.digest, 'accepted_var_name']
                            #     old_docker_var_name = task.runtime['docker'].expr.referee.name
                            #     default_docker_image = str(var.expr)
                            #     self.convert_task_input_to_no_default(document, task, old_docker_var_name, docker_var_name)
                            #     self.change_task_runtime_docker_var_name(document, task, docker_var_name)
                            #     self.change_input_and_update_callers(tree_parser, task, old_docker_var_name, docker_var_name, default_docker_image)
                            # else:
                            #     docker_var_name = task.runtime['docker'].expr.referee.name
                            #     self.convert_task_input_to_no_default(document, task, docker_var_name, docker_var_name)
                        else:
                            assert False, f"Unknown type of expression: {type(var.expr)}, {task}"
                    else:
                        assert False, f"Variable is not a decl: {str(var)}, {task}"
                case _:
                    assert False, f"Unknown type of docker image variable: {type(task.runtime['docker'])}, {task}"
        return

    def get_document_meta(self, document):
        if type(document) != WDL.Tree.Document:
            # Walk up the parents until we find the WDL.Tree.Document:
            parent = document
            while parent.parent and type(parent) != WDL.Tree.Document:
                parent = parent.parent
            assert type(parent) == WDL.Tree.Document, f"Could not find WDL.Tree.Document in parents of {document}"
            document = parent

        if document.pos.abspath in self.modified_docs:
            return self.modified_docs[document.pos.abspath]
        else:
            assert False, f"Document not found in modified_docs: {document.pos.abspath}"

    def get_new_docker_var_name(self, task):
        row = self.var_name_df[self.var_name_df.index == task.digest]
        assert row.shape[0] == 1, f"Expected to find a single row for task: {task.name} [{task.digest}]"
        old_var_name = row['var_name'][0]
        new_var_name = row['accepted_var_name'][0]
        if row['var_name_override'][0]:
            assert old_var_name, f"For a variable override there should be an old variable name: {task.name} [{task.digest}], {new_var_name}"
            log.info(f"Changing docker image variable name, from '{old_var_name}' to '{new_var_name}', for task: '{task.name}' ['{task.digest}']")
        # assert old_var_name, f"Expected to find an old variable name for task: {task.name} [{task.digest}], {new_var_name}"
        return old_var_name, new_var_name

    def get_l_justify_inputs(self, meta, task):
        """Get the minimum number of spaces on input lines to left justify the input variables.
        ' ' * (len(line) - len(line.lstrip())) does not work because sometimes the last line is weirdly indented."""
        match(type(task)):
            case WDL.Tree.Task | WDL.Tree.Workflow:
                lines = [x.pos.line for x in task.inputs]
            case WDL.Tree.Call:
                lines = [task.inputs[x].pos.line for x in task.inputs]
            case _:
                assert False, f"Unknown type of task: {type(task)}"
        return ' ' * min([len(meta.original_lines[x-1]) - len(meta.original_lines[x-1].lstrip()) for x in lines])

    def patch_task_input(self, task, old_var_name, new_var_name, keep_default_val=False, var_type="String"):
        """Patch task inputs to replace/rename or add a new variable. 
        input {
            String old_var_name = "broadinstitute/gatk"    =>    String new_var_name
            String docker = "broadinstitute/gatk"          =>    String gatk_docker
            String docker = "broadinstitute/gatk"          =>    String gatk_docker = "broadinstitute/gatk"
            String docker = "broadinstitute/gatk"          =>    String docker

        """
        meta = self.get_document_meta(task)
        assert type(task) == WDL.Tree.Task or type(task) == WDL.Tree.Workflow, f"Expected task type to be WDL.Tree.Task or Workflow, found: {type(task)}"
        default_val = ''
        var_line_pos, var_obj = self.patch_get_input_or_last_input_line(task, old_var_name)
        if keep_default_val and old_var_name and var_obj and var_obj.expr:
            default_val = ' = ' + str(var_obj.expr)
            assert var_type in str(type(var_obj.type)), f"Variable type does not match expected type: {var_type}, {var_obj.type}"
        # Patch the input line (or add a new one):
        line = meta.original_lines[var_line_pos - 1]
        l_justify = self.get_l_justify_inputs(meta, task)
        # new_line = f"{l_justify}{var_type} {new_var_name}{default_val} [task_input: {type(task)}]\n"
        new_line = f"{l_justify}{var_type} {new_var_name}{default_val}\n"
        if var_obj:
            # Replace the old line with the new line:
            meta.modified_lines[var_line_pos - 1] = new_line
        else:
            # Append the new line to the end of the inputs:
            if new_line not in meta.modified_lines[var_line_pos - 1]:
                meta.modified_lines[var_line_pos - 1 ] = meta.modified_lines[var_line_pos - 1 ] + new_line
            pass # TODO: Finish with a mechanism to avoid repeats

    def patch_task_runtime_parameter(self, task, runtime_param_name, new_value, keep_old_value=False):
        """Patch runtime_param_name to have a new variable (new_value) or create a new variable if runtime_param_name is not found
        keep_default_val=True (typically not used) will keep the default value whether it is a String or a variable reference
        runtime {
            runtime_param_name: gatk_docker    ==>    runtime_param_name: new_value
            docker: gatk_docker                ==>    docker: new_value
            docker: "us.gcr.io/ubuntu:1.1.2"   ==>    docker: new_value
        """
        meta = self.get_document_meta(task)
        assert type(task) == WDL.Tree.Task, f"Expected task type to be WDL.Tree.Task, found: {type(task)}"
        default_val = ': ' + new_value
        var_obj = None
        if runtime_param_name and runtime_param_name in task.runtime:
            var_obj = task.runtime[runtime_param_name]
            var_line_pos = var_obj.pos.line
            if keep_old_value and type(var_obj) == WDL.Expr.String:
                default_val = ': ' + str(var_obj)
            elif keep_old_value and type(var_obj) == WDL.Expr.Get: 
                default_val = ': ' + str(var_obj.expr.referee.name)
        else:
            var_line_pos = max([task.runtime[x].pos.end_line for x in task.runtime])
        # Patch the input line (or add a new one):
        line = meta.original_lines[var_line_pos - 1]
        l_justify = self.get_l_justify_inputs(meta, task)
        # new_line = f"{l_justify}{runtime_param_name}{default_val} [task_runtime]\n"
        new_line = f"{l_justify}{runtime_param_name}{default_val}\n"
        assert runtime_param_name != None, f"Expected runtime_param_name to be defined"
        if var_obj:
            # Replace the old line with the new line:
            meta.modified_lines[var_line_pos - 1] = new_line
        else:
            # Append the new line to the end of the runtime inputs:
            # TODO: needs better fix
            if new_line not in meta.modified_lines[var_line_pos - 1]:
                meta.modified_lines[var_line_pos - 1 ] = meta.modified_lines[var_line_pos - 1 ] + new_line
            pass

    def patch_caller_input(self, task, old_input_name, new_input_name, new_input_value, keep_default_val=False, var_type="String"):
        """Patch caller inputs to replace/rename or add a new variable, and set a new value which should be from the task's inputs/namespace
        call Workflow.Task as Alias {
            input:
                old_input_name = default_val    =>    new_input_name = new_input_value
                gatk_docker = gatk_docker       =>    gatk_docker = references.gatk_docker
                docker = references.gatk_docker =>    gatk_docker = gatk_docker
        
        Special cases for WDL.Tree.Call - inputs to a call are separated by a comma:
            input_vcf = input_vcf     =>    input_vcf = input_vcf,    [must add a comma to an previously last input line]
                                            gatk_docker = gatk_docker,
                                            python_docker = python_docker
        
            input_vcf = input_vcf,    =>    input_vcf = input_vcf,
            docker = docker,                gatk_docker = gatk_docker,  [must add a comma on a replacement when not the last input line]
            input_ref = input_ref,                
        """
        meta = self.get_document_meta(task)
        assert type(task) == WDL.Tree.Call, f"Expected task type to be WDL.Tree.Call, found: {type(task)}"
        default_val = ' = ' + new_input_value
        var_line_pos, var_obj = self.patch_get_input_or_last_input_line(task, old_input_name)
        if keep_default_val and old_input_name and var_obj and var_obj.expr:
            default_val = ' = ' + str(var_obj.expr)
            assert var_type in str(type(var_obj.type)), f"Variable type does not match expected type: {var_type}, {var_obj.type}"
        # Patch the input line (or add a new one):
        line = meta.original_lines[var_line_pos - 1]
        l_justify = self.get_l_justify_inputs(meta, task)
        # new_line = f"{l_justify}{new_input_name}{default_val} [caller input]\n"
        new_line = f"{l_justify}{new_input_name}{default_val}\n"
        if var_obj:
            if var_line_pos <= max([task.inputs[x].pos.line for x in task.inputs]):
                # This is not the last line, so add a comma:
                meta.modified_lines[var_line_pos - 1] = new_line.rstrip() + ',\n'
            else:
                meta.modified_lines[var_line_pos - 1] = new_line
        else:
            # Confirm that the line we're appending to has a ,\n at the end:
            if meta.modified_lines[var_line_pos - 1].rstrip()[-1] != ',':
                meta.modified_lines[var_line_pos - 1] = meta.modified_lines[var_line_pos - 1].rstrip() + ',\n'
            # Only append if the new_line isn't already there:
            if new_line not in meta.modified_lines[var_line_pos - 1] and (new_line.rstrip() + ',\n') not in meta.modified_lines[var_line_pos - 1]:
                meta.modified_lines[var_line_pos - 1 ] = meta.modified_lines[var_line_pos - 1 ] + new_line

    def patch_get_input_or_last_input_line(self, task, old_var_name):
        input_line_pos = None
        input_var = None
        last_input_line = None
        if type(task) == WDL.Tree.Task or type(task) == WDL.Tree.Workflow:
            matching_input = [x for x in task.inputs if x.name == old_var_name]
            if len(matching_input) >= 1:  # Matching input found
                assert len(matching_input) == 1, f"Found multiple matching inputs: {matching_input}"
                assert matching_input[0].pos.line == matching_input[0].pos.end_line, f"Found input with multiple lines: {matching_input[0]}"
                input_line_pos = matching_input[0].pos.line
                input_var = matching_input[0]
            else:  # No matching input, add a new line at the end of the inputs
                last_input_line = max([x.pos.end_line for x in task.inputs])
        elif type(task) == WDL.Tree.Call:
            matching_input = [v for k,v, in task.inputs.items() if k == old_var_name]
            if len(matching_input) >= 1:  # Matching input found
                assert len(matching_input) == 1, f"Found multiple matching inputs: {matching_input}"
                assert matching_input[0].pos.line == matching_input[0].pos.end_line, f"Found input with multiple lines: {matching_input[0]}"
                input_line_pos = matching_input[0].pos.line
                input_var = matching_input[0]
            else:  # No matching input, add a new line at the end of the inputs
                last_input_line = max([v.pos.end_line for k,v, in task.inputs.items()])
        if input_line_pos:
            return input_line_pos, input_var
        elif last_input_line:
            return last_input_line, None
        else:
            assert False, f"Could not find input line for {task.name} {old_var_name}"

    def patch_update_caller_hierarchy(self, tree_parser, patched_task, old_docker_var_name, new_docker_var_name):
        """Tasks inside WDLs are called by other tasks (first with a WDL.Tree.Call, caller) and then with a call hierarchy up various task/sub-workflow levels
        Each caller, task, and workflow needs to pass the appropriate input variables to this patched_task. Or the caller needs to use a fully-qualified variable to set the input variable
        A patched_task may have many callers (represented in the tree_parser.tree as a [caller_digest]__[library_call (our patched_task.digest)])
        For cases where the docker variable is being changed the old_docker_var_name != new_docker_var_name (changing generic 'docker' to specific 'gatk_docker')
        For cases where the docker variable wasn't in the caller old_docker_var_name can either be set or be None, but new_docker_var_name needs to be set (adding inputs that were not passed
        i.e., cases where the caller set a default instead of taking an input from the user)
        """
        callers = [x for x in self.find_task_calls(patched_task, tree_parser)]
        for direct_caller in callers:
            task = tree_parser.tree_map[direct_caller['digest']]
            assert type(task) == WDL.Tree.Call, f"Expected first task type to be WDL.Tree.Call, found: {type(task)}"

            # Now we can update all tasks/workflows that contain this call
            for call_parent_digest in [direct_caller['digest']] + direct_caller['parent_digest']:
            # for call_parent_digest in direct_caller['parent_digest']:
                assert call_parent_digest in tree_parser.tree_map, f"Could not find call parent {call_parent_digest} in tree_map"
                parent = tree_parser.tree_map[call_parent_digest]
                parent_node = [x for x in tree_parser.tree if x['digest'] == call_parent_digest]
                assert len(parent_node) == 1, f"Expected to find a single node with digest {call_parent_digest}, found {len(parent_node)}"
                parent_node = parent_node[0]
                if parent_node['type'] == 'scatter':
                    continue  # Scatters are just logical entities they shouldn't modify inputs
                
                match(type(parent)):
                    case WDL.Tree.Task | WDL.Tree.Workflow:
                        self.patch_task_input(parent, old_docker_var_name, new_docker_var_name)
                    case WDL.Tree.Call:
                        # The caller is special because it has a section to format the input parameters we need to modify:
                        self.patch_caller_input(parent, old_docker_var_name, new_docker_var_name, new_docker_var_name)
                        # Callers have a digest that is TaskOrWorkflowDigest__CallerDigest, so we need to find the workflow containing this call:
                        task_or_workflow_digest = call_parent_digest.split('__')[0]
                        caller_parent_obj = tree_parser.caller_tree_map[task_or_workflow_digest]
                        # search_parent = parent.parent
                        # while search_parent.parent:
                        #     if hasattr(search_parent, 'digest') and search_parent.digest == workflow_digest:
                        #         break
                        #     search_parent = search_parent.parent
                        assert caller_parent_obj.digest == task_or_workflow_digest, f"Could not find task or workflow with digest {task_or_workflow_digest} in tree_map"
                        self.patch_task_input(caller_parent_obj, old_docker_var_name, new_docker_var_name)

                        # self.patch_task_input(parent, old_docker_var_name, new_docker_var_name)
        return

    def store_docker_default(self, document, task, docker_variable_name, default_docker_image):
        uniq_task_name = f"{document.pos.abspath}_{task.name}_{task.pos.line}"
        # TODO: Check if this variable already exists with another name, if so use that name
        self.modified_tasks[uniq_task_name] = [docker_variable_name, default_docker_image]
        return uniq_task_name

    def find_task_calls(self, task, tree_parser):
        for call in tree_parser.tree:
            if call['type'] != 'task':
                continue
            # Digest sanity check:
            if call['callee_name'] == task.name and call['callee_pos'].abspath == task.pos.abspath:
                if 'callee_digest' in call and call['callee_digest'] != task.digest:
                    assert False, f"Digest mismatch: {call['callee_digest']} != {task.digest}"
            if 'callee_digest' in call and call['callee_digest'] == task.digest:
                yield call
