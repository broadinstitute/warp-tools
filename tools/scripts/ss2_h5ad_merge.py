#!/usr/bin/env python3

import argparse
import anndata as ad
import os

def merge_h5ads(input_h5ad_files, output_h5ad_file):
    """
    Merge multiple H5AD files into a single H5AD file.
    
    Args:
        input_h5ad_files (list of str): Paths to input H5AD files.
        output_h5ad_file (str): Path to output merged H5AD file.
    """
    # Initialize an empty list to store the loaded AnnData objects
    adata_list = []

    # Load each H5AD file and append it to the list
    for h5ad_file in input_h5ad_files:
        print(f"Loading {h5ad_file}...")
        adata = ad.read_h5ad(h5ad_file)
        adata_list.append(adata)
    
    # Concatenate all the AnnData objects along the observation (cell) axis
    print("Merging H5AD files...")
    merged_adata = ad.concat(adata_list, axis=0)
    
    # Save the merged AnnData object to a new H5AD file
    print(f"Saving merged H5AD file to {output_h5ad_file}...")
    merged_adata.write_h5ad(output_h5ad_file)
    print("Merge complete!")

def main():
    description = "Merge multiple H5AD files into a single H5AD file."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-h5ad-files',
                        dest='input_h5ad_files',
                        nargs="+",
                        required=True,
                        help="Paths to input H5AD files")
    parser.add_argument('--output-h5ad-file',
                        dest='output_h5ad_file',
                        required=True,
                        help="Path to output merged H5AD file")
    args = parser.parse_args()

    merge_h5ads(args.input_h5ad_files, args.output_h5ad_file)

if __name__ == '__main__':
    main()

