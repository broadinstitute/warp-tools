import argparse
import anndata as ad
import pandas as pd

def merge_h5ad_files(args):
    """This function merges multiple H5AD files into one, preserving metadata."""
    h5ad_files = args.input_h5ad_files
    adatas = [ad.read_h5ad(f) for f in h5ad_files]

    # Concatenate the AnnData objects
    merged_adata = ad.concat(adatas, axis=0)

    # Add global attributes (uns)
    attrDict = dict()
    attrDict['batch_id'] = args.batch_id
    attrDict['pipeline_version'] = args.pipeline_version
    if args.batch_name is not None:
        attrDict['batch_name'] = args.batch_name
    if args.library is not None:
        attrDict['library_preparation_protocol.library_construction_approach'] = args.library
    if args.species is not None:
        attrDict['donor_organism.genus_species'] = args.species
    if args.organ is not None:
        attrDict['specimen_from_organism.organ'] = args.organ
    if args.project_id is not None:
        attrDict['project.provenance.document_id'] = args.project_id
    if args.project_name is not None:
        attrDict['project.project_core.project_short_name'] = args.project_name
    
    # Assign these attributes to the merged AnnData object
    merged_adata.uns.update(attrDict)
    
    # Save the merged H5AD file
    merged_adata.write_h5ad(args.output_h5ad_file)

def main():
    description = """Merge the outputs of multiple SS2 pipeline runs into a single H5AD file"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-h5ad-files',
                        dest='input_h5ad_files',
                        nargs="+",
                        required=True,
                        help="Path to input H5AD files")
    parser.add_argument('--output-h5ad-file',
                        dest='output_h5ad_file',
                        required=True,
                        help="Path to output H5AD file")
    parser.add_argument('--batch_id',
                        dest='batch_id',
                        required=True,
                        help="Batch id for output H5AD")
    parser.add_argument('--batch_name',
                        dest='batch_name',
                        help='User provided batch name for output H5AD')
    parser.add_argument('--project_id',
                        dest='project_id',
                        help='User provided project id for output H5AD')
    parser.add_argument('--project_name',
                        dest='project_name',
                        help='User provided project name for output H5AD')
    parser.add_argument('--library',
                        dest='library',
                        help='User provided library for output H5AD')
    parser.add_argument('--species',
                        dest='species',
                        help='User provided species for output H5AD')
    parser.add_argument('--organ',
                        dest='organ',
                        help='User provided organ for output H5AD')
    parser.add_argument('--pipeline_version',
                        dest='pipeline_version',
                        required=True,
                        help='Multisample SS2 version')
    args = parser.parse_args()

    merge_h5ad_files(args)

if __name__ == '__main__':
    main()

