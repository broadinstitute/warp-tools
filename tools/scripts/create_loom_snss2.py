import argparse
import csv
import numpy as np
import pandas as pd
import anndata
import loompy
from scipy.sparse import csc_matrix

def generate_col_attr(args):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to loom file
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
                numeric_value = float(value.strip("%")) / 100
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

    for i in range(numeric_field_names.shape[0]):
        name = numeric_field_names[i]
        data = np.array([sorted_numeric_values])[:,i]
        col_attrs[name] = data
    string_field_names = np.array(sorted_string_labels)
    for i in range(string_field_names.shape[0]):
        name = string_field_names[i]
        data = np.array([sorted_string_values])[:,i]
        col_attrs[name] = data

    return col_attrs


def generate_csr_sparse_coo(expr):
    expr_coo = csc_matrix(expr[:])
    return expr_coo

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

def create_h5ad_file(args):
    col_attrs = generate_col_attr(args)
    row_attrs, intron_expr_csr_coo, exon_expr_csr_coo = generate_row_attr_and_matrix(args.count_results_file)

    col_attrs_df = pd.DataFrame(col_attrs, index=[args.input_id])
    row_attrs_df = pd.DataFrame(row_attrs, index=row_attrs['Gene'])

    adata = anndata.AnnData(X=exon_expr_csr_coo, obs=col_attrs_df, var=row_attrs_df)
    adata.layers['intron_counts'] = intron_expr_csr_coo.tocsc()
    adata.write(args.output_h5ad_path)

def compare_loom_to_h5ad(loom_path, h5ad_path):
    loom_file = loompy.connect(loom_path)
    adata = anndata.read_h5ad(h5ad_path)

    loom_data = loom_file[:, :]
    h5ad_data = adata.X.toarray()

    assert np.array_equal(loom_data, h5ad_data), "Data in LOOM and H5AD files differ!"

    loom_file.close()
    print("LOOM and H5AD files are equivalent.")

def main():
    parser = argparse.ArgumentParser(description="Convert pipeline outputs to H5AD and compare with LOOM.")
    parser.add_argument('--qc_files', nargs="+", help='QC files')
    parser.add_argument('--count_results', dest="count_results_file", help='Path to count results')
    parser.add_argument('--output_h5ad_path', help='Path to save the H5AD file')
    parser.add_argument('--input_id', help='Sample or cell ID')
    parser.add_argument('--pipeline_version', default="Unknown sample", help='Pipeline version')
    parser.add_argument('--compare_to_loom', help='Path to the existing LOOM file for comparison')

    args = parser.parse_args()

    create_h5ad_file(args)

    if args.compare_to_loom:
        compare_loom_to_h5ad(args.compare_to_loom, args.output_h5ad_path)

if __name__ == '__main__':
    main()
