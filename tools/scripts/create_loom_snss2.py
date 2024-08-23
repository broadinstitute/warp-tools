import numpy as np
import scipy.sparse as sc
import pandas as pd
import anndata as ad
import csv

def generate_col_attr(args):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to loom file
    Args:
        args (Namespace): Argument namespace containing paths and IDs.
    Returns:
        col_attrs (dict): Column attributes dictionary.
    """
    # read the QC values
    qc_path = [p for p in args.qc_files if p.endswith(".csv")][0]
    qc_values = pd.read_csv(qc_path).values
    n_metrics_files = len(qc_values) - 2
    metadata_labels = qc_values[0, 1:]
    metadata_values = []
    for index_val in range(len(metadata_labels)):
        arr_metric_val = []
        for index_file in range(n_metrics_files):
            arr_metric_val.append(qc_values[index_file + 2, index_val + 1])
        metadata_values.append(','.join(s for s in arr_metric_val if s))

    cell_id = args.input_id
    string_metadata = {}
    numeric_metadata = {}

    for label, value in zip(metadata_labels, metadata_values):
        # Check if value is numeric or string
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
    sorted_string_labels = sorted(string_metadata.keys())
    sorted_string_values = [string_metadata[m] for m in sorted_string_labels]
    sorted_numeric_labels = sorted(numeric_metadata.keys())
    sorted_numeric_values = [numeric_metadata[m] for m in sorted_numeric_labels]

    # Column attributes
    col_attrs = {}
    col_attrs["cell_names"] = [cell_id]
    col_attrs["CellID"] = [cell_id]
    col_attrs['input_id'] = [args.input_id]

    numeric_field_names = np.array(sorted_numeric_labels)
    for i in range(numeric_field_names.shape[0]):
        name = numeric_field_names[i]
        data = np.array([sorted_numeric_values])[:, i]
        col_attrs[name] = data

    string_field_names = np.array(sorted_string_labels)
    for i in range(string_field_names.shape[0]):
        name = string_field_names[i]
        data = np.array([sorted_string_values])[:, i]
        col_attrs[name] = data

    return col_attrs

def generate_csr_sparse_coo(expr):
    """Converts a given expression matrix to CSR format if it is in COO format.
    Args:
        expr (scipy.sparse.coo_matrix or np.ndarray): Expression matrix.
    Returns:
        scipy.sparse.csr_matrix: CSR format matrix.
    """
    if isinstance(expr, np.ndarray):
        expr = sc.coo_matrix(expr)
    if not isinstance(expr, sc.coo_matrix):
        raise ValueError("Input must be a COO matrix or an ndarray.")
    return expr.tocsr()

def create_h5ad_file(args):
    """Creates an H5AD file from the provided QC and count results files.
    Args:
        args (Namespace): Argument namespace containing paths and IDs.
    """
    # Process QC and count results
    col_attrs = generate_col_attr(args)

    # Read count results
    count_results = pd.read_csv(args.count_results, delimiter='\t')
    gene_ids = count_results['gene_id'].tolist()
    gene_names = count_results['gene_name'].tolist()
    intron_counts = count_results['introns'].tolist()
    exon_counts = count_results['exons'].tolist()
    intron_lengths = count_results['intron_length'].tolist()
    exon_lengths = count_results['exon_length'].tolist()

    # Convert counts to CSR format
    intron_expression_csr_coo = generate_csr_sparse_coo([intron_counts])
    exon_expression_csr_coo = generate_csr_sparse_coo([exon_counts])

    # Create AnnData object
    adata = ad.AnnData(
        X=intron_expression_csr_coo,
        obs=pd.DataFrame(index=gene_ids),
        var=pd.DataFrame(index=gene_names),
        uns={'pipeline_version': args.pipeline_version}
    )

    # Adding column attributes
    for key, value in col_attrs.items():
        adata.obs[key] = value

    # Convert matrix to CSR format if needed
    if isinstance(adata.X, sc.coo_matrix):
        adata.X = adata.X.tocsr()

    # Write AnnData object to H5AD file
    adata.write(args.output_h5ad_path)

