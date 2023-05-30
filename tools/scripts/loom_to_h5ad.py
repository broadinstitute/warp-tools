import argparse
import anndata as ad
import loompy
import os

def convertloomtoh5ad(args):
    """ Converts loom file format to h5ad. 

    Args:
        loom_path (str): path to loom file.  
        output_path (str): output path directory to save h5ad file in. 
    
    """

    file_name = os.path.basename(args.loom_file).split(".")[0] 

    #reads loom file as sparse csr matrix
    adata = ad.read_loom(args.loom_file, sparse=True)
    
    #connect to loom file to save the global attributes in the h5ad file
    ds = loompy.connect(args.loom_file)
    
    #loop through global attributes and save in adata
    for global_attr, value in ds.attrs.items():
        adata.uns[global_attr] = value
    ds.close()
    
    adata.write(os.path.join(args.output_path, file_name + ".h5ad"))

def main():
    description = """This script converts a loom file (http://linnarssonlab.org/loompy/index.html) into h5ad (https://anndata.readthedocs.io/en/latest/) format."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument(
        "--loom_file",
        dest="loom_file",
        required=True,
        help="A loom file is needed as input.",
    )

    parser.add_argument(
        "--output_path",
        dest="output_path",
        required=True,
        help="Output path to save the h5ad file in.",
    )

    args = parser.parse_args()
    convertloomtoh5ad(args)   

if __name__ == "__main__":
    main()
