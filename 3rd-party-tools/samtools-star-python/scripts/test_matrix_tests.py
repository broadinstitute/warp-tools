import pytest
import pandas as pd
#import sys
#sys.path.append('/warp-tools/3rd-party-tools/star-merge-npz/scripts/combine_shard_metrics.py')
import os

from combine_shard_metrics import merge_matrices

@pytest.fixture
def example_input_files():
    # Get the absolute path of the current script
    absolute_script_path = os.path.abspath(__file__)
    print("The absolute path is :")
    print(absolute_script_path)

    # Get the parent directory of the script
    parent_directory = os.path.dirname(absolute_script_path)
    parent_directory_absolute = os.path.abspath(parent_directory)
    print("Absolute directory is: ")
    print(parent_directory_absolute)

    return {
        "summary_file": parent_directory_absolute + "/input/Human_v2_snRNA_summary.txt",
        "align_file": parent_directory_absolute + "/input/Human_v2_snRNA_align_features.txt",
        "cell_reads": parent_directory_absolute + "/input/Human_v2_snRNA_cell_reads.txt",
        "counting_mode": "sn_rna",
        "uniform_barcodes": parent_directory_absolute + "/input/barcodes.tsv",
        "uniform_mtx": parent_directory_absolute + "/input/matrix.mtx",
        "expected_cells": "3000"
    }

def test_merge_matrices_column_existence(example_input_files):
    df = merge_matrices(**example_input_files)

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