import pytest
import pandas as pd
#import sys
#sys.path.append('/warp-tools/3rd-party-tools/star-merge-npz/scripts/combine_shard_metrics.py')

from .. import combine_shard_metrics

@pytest.fixture
def example_input_files():
    return {
        "summary_file": "input/Human_v2_snRNA_summary.txt",
        "align_file": "input/Human_v2_snRNA_align_feautres.txt",
        "cell_reads": "input/Human_v2_snRNA_cell_reads.txt",
        "counting_mode": "sn_rna",
        "base_name": "output",
        "uniform_barcodes": "input/barcodes.tsv",
        "uniform_mtx": "input.matrix.mtx"
    }

def test_merge_matrices_column_existence(example_input_files):
    df = merge_matrices(example_input_files)

    # Define expected columns
    expected_columns = [
        "number_of_reads",
        "sequencing_saturation",
        "fraction_of_unique_reads_mapped_to_genome",
        "fraction_of_unique_and_multiple_reads_mapped_to_genome",
        "fraction_of_reads_with_Q30_bases_in_rna",
        "fraction_of_reads_with_Q30_bases_in_cb_and_umi",
        "fraction_of_reads_with_valid_barcodes",
        "reads_mapped_antisense_to_gene",
        "reads_mapped_confidently_exonic",
        "reads_mapped_confidently_to_genome",
        "reads_mapped_confidently_to_intronic_regions",
        "reads_mapped_confidently_to_transcriptome",
        "estimated_cells",
        "umis_in_cells",
        "mean_umi_per_cell",
        "median_umi_per_cell",
        "unique_reads_in_cells_mapped_to_gene",
        "fraction_of_unique_reads_in_cells",
        "mean_reads_per_cell",
        "median_reads_per_cell",
        "mean_gene_per_cell",
        "median_gene_per_cell",
        "total_genes_unique_detected"
    ]

    # Ensure all expected columns are present in the dataframe
    for column in expected_columns:
        assert column in df.columns, f"Column '{column}' is missing from the dataframe"