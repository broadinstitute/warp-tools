import argparse
import csv
import os
import numpy as np
import scipy as sc
import anndata as ad
import pandas as pd

def generate_col_attr(args):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to H5AD file
    Args:
        qc_path (str): path to the QCs csv
    Retruns:
         col_attrs: dict
    """
    # read the QC values
    qc_path = [p for p in args.qc_files if p.endswith(".csv")][0]
    with open(qc_path, 'r') as f:
        qc_values = [row for row in csv.reader(f)]
        n_metrics_files = len(qc_values)-2
        metadata_labels = qc_values[0][1:]
        metadata_values = []
        for index_val in range(len(metadata_labels)):
          arr_metric_val = []
          for index_file in range(n_metrics_files):
            arr_metric_val.append(qc_values[index_file+2][index_val+1])
          metadata_values.append(','.join(s for s in arr_metric_val if s))

    cell_id = args.input_id

    string_metadata = {}
    numeric_metadata = {}

    for label, value in zip(metadata_labels, metadata_values):

        # See if this is numeric or string
        numeric_value = None
        try:
            numeric_value = float(value)
        except ValueError:
            try:
                numeric_value = float(value.strip("%"))/100
            except ValueError:
                pass

        if numeric_value is not None:
            numeric_metadata[label] = numeric_value
        else:
            string_metadata[label] = value

    # Metrics
    # Write the string and numeric metadata separately
    sorted_string_labels = sorted(string_metadata.keys())
    sorted_string_values = [string_metadata[m] for m in sorted_string_labels]
    sorted_numeric_labels = sorted(numeric_metadata.keys())
    sorted_numeric_values = [numeric_metadata[m] for m in sorted_numeric_labels]

    # Column attributes
    col_attrs = dict()
    col_attrs["cell_names"] = [cell_id]
    col_attrs["CellID"] = [cell_id]
    col_attrs['input_id'] = [args.input_id]

    numeric_field_names = np.array(sorted_numeric_labels[:])

    for i in range(0, numeric_field_names.shape[0]):
        name = numeric_field_names[i]
        data = np.array([sorted_numeric_values])[:,i]
        col_attrs[name] = data
    string_field_names = np.array(sorted_string_labels)
    for i in range(0, string_field_names.shape[0]):
        name = string_field_names[i]
        data = np.array([sorted_string_values])[:,i]
        col_attrs[name] = data

    return col_attrs


def generate_csr_sparse_coo(expr):

    nrows, ncols = np.shape(expr)
    expr_coo = sc.sparse.coo_matrix(expr[:])
    xcoord = []
    ycoord = []
    value = []

    for k in range(0, expr_coo.data.shape[0]):
        xcoord.append(expr_coo.row[k])
        ycoord.append(expr_coo.col[k])
        value.append(expr_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(ncols, nrows))
    return expr_sp_t


def generate_row_attr_and_matrix(count_results_path):
    """Converts the Smart Seq2 intron and exon counts file to
        csr_coos and row attribute dictionary
    Args:
        input_path (str): file where the SS2 pipeline expression counts are

    Return:
        row_attrs: dict of attributes
        intron_expression_csr_coo: crs coo matrix for introns
        exon_expression_csr_coo: crs coo matrix for exons
    """
    reader = csv.DictReader(open(count_results_path), delimiter="\t")

    intron_expression_values = {}
    exon_expression_values = {}
    intron_lengths = {}
    exon_lengths = {}

    gene_ids = []
    gene_names = []
    count_values = {}
    for row in reader:
        intron_expression_values[row["gene_id"]] = int(row["introns"])
        exon_expression_values[row["gene_id"]] = int(row["exons"])
        intron_lengths[row["gene_id"]] = int(row["intron_length"])
        exon_lengths[row["gene_id"]] = int(row["exon_length"])
        gene_ids.append(row['gene_id'])
        gene_names.append(row['gene_name'])

    intron_counts = [intron_expression_values[g] for g in gene_ids]
    exon_counts  = [exon_expression_values[g] for g in gene_ids]
    intron_lengths = [intron_lengths[g] for g in gene_ids]
    exon_lengths = [exon_lengths[g] for g in gene_ids]

    intron_expression_csr_coo = generate_csr_sparse_coo([intron_counts])
    exon_expression_csr_coo = generate_csr_sparse_coo([exon_counts])

    row_attrs = { "ensembl_ids"  : np.array(gene_ids),
                 "gene_names"    : np.array(gene_names),
                 "Gene"          : np.array(gene_ids),
                 "intron_lengths": np.array(intron_lengths),
                 "exon_lengths"  : np.array(exon_lengths)
                }

    return row_attrs, intron_expression_csr_coo, exon_expression_csr_coo


def create_h5ad_files(args):
    """This function creates the H5AD file or folder structure in output_h5ad_path in
       format file_format, with input_id from the input folder analysis_output_path
    Args:
        input_id (str): sample or cell id
        qc_analysis_output_files_string (str): a string with the file names in the QCGroup of SS2
            pipeline output, separated by commas
        rsem_genes_results_file (str): the file for the expression count
        output_h5ad_path (str): location of the output H5AD file
    """
    # generate a dictionary of column attributes
    col_attrs =  generate_col_attr(args)

    # add the expression count matrix data
    # generate a dictionary of row attributes
    row_attrs, intron_expr_csr_coo, exon_expr_csr_coo = generate_row_attr_and_matrix(args.count_results_file)

    # Convert the matrices to CSR format
    exon_counts = exon_expr_csr_coo.transpose().tocsr()
    intron_counts = intron_expr_csr_coo.transpose().tocsr()

    gene_attrs = row_attrs
    cell_attrs = col_attrs
    adata = ad.AnnData(exon_counts)
    attrDict = dict()
    attrDict['pipeline_version'] = args.pipeline_version

    # importing the cell and gene metrics as dataframes
    dc = pd.DataFrame(cell_attrs)
    dg = pd.DataFrame(gene_attrs)

    # assign the cell and gene attrs to the obs and var field respectively
    adata.obs = dc
    adata.var = dg

    # assign global attrs to unstructured data
    adata.uns = attrDict
    adata.obs_names = [x for x in adata.obs['CellID']]

    # set variable names
    adata.var_names = [x for x in adata.var['Gene']]

    # set the layers to the intron counts
    adata.layers['intron_counts'] = intron_counts

    # Write the h5ad file
    adata.write_h5ad(args.output_h5ad_path + ".h5ad")

def main():
    description = """This script converts the some of the SmartSeq2 pipeline outputs in to
                   H5AD format. This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--qc_files',
                        dest="qc_files",
                        nargs = "+",
                        help=('the grouped QC files from the GroupQCOutputs task of SS2 '
                              'Single Sample workflow'))

    parser.add_argument('--count_results',
                        dest="count_results_file",
                        help='path to the intronic and exonic count and FPKM to be added to the H5AD')

    parser.add_argument('--output_h5ad_path',
                        dest="output_h5ad_path",
                        help='path where the h5ad file is to be created')

    parser.add_argument('--input_id',
                        dest="input_id",
                        help='the sample name in the bundle')

    parser.add_argument("--input_name",
                        dest="input_name",
                        help= "sequencing_input.biomaterial_core.biomaterial_id in HCA metadata, defined by the user",
    )

    parser.add_argument("--input_id_metadata_field",
                        dest="input_id_metadata_field",
                        help= "sequencing_process.provenance.document_id: [UUID] defined by the user",
    )

    parser.add_argument("--input_name_metadata_field",
                        dest="input_name_metadata_field",
                        help= "sequencing_input.biomaterial_core.biomaterial_id defined by the user",
    )

    parser.add_argument('--pipeline_version',
                        default="Unknown sample",
                        help='the version of SS2 used to generate data')

    args = parser.parse_args()

    create_h5ad_files(args)


if __name__ == '__main__':
    main()

