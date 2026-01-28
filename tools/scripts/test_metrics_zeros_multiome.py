#!/usr/bin/env python3
"""
Script for checking zero metrics in Optimus output files.

This script checks for columns/metrics that are all zeros in gene metrics,
cell metrics, and library metrics files. Tests are skipped if the corresponding
file argument is not provided.

Currently this only checks for zeros expected in Multiome mode, so if this is run on Optimus outputs it may show errors for expected zero metrics.

"""

import argparse
import pandas as pd
import sys
import os

# Metrics that are expected to be all zeros (will verify they ARE zero)
GENE_METRICS_EXPECTED_ZERO = [
    "reads_mapped_exonic_as",
    "noise_reads",
    "reads_mapped_intronic_as",
    "duplicate_reads",
    "antisense_reads",
]

CELL_METRICS_EXPECTED_ZERO = [
    "noise_reads",
    "duplicate_reads",
    "antisense_reads",
    "reads_unmapped",
    "reads_mapped_too_many_loci",
]

LIBRARY_METRICS_EXPECTED_ZERO = [
]

LIBRARY_METRICS_EXPECTED_ZERO_SCRNA_MODE = [
    "reads_mapped_antisense_to_gene",
    "reads_mapped_confidently_to_intronic_regions",
    "percent_intronic_reads",
]

def _load_csv(file_path: str) -> pd.DataFrame:
    """Load a CSV or TSV file, auto-detecting the delimiter."""
    try:
        return pd.read_csv(file_path, sep=None, engine='python')
    except Exception:
        return pd.read_csv(file_path)


def _find_all_zero_columns(df: pd.DataFrame, expected_zero: list) -> tuple:
    """
    Find numeric columns where all values are zero.
    Returns (unexpected_zero_cols, expected_but_not_zero_cols, expected_and_zero_cols).
    """
    unexpected_zero = []
    expected_but_not_zero = []
    expected_and_zero = []
    
    numeric_cols = df.select_dtypes(include='number').columns
    
    for col in numeric_cols:
        is_all_zero = (df[col] == 0).all()
        
        if col in expected_zero:
            if is_all_zero:
                expected_and_zero.append(col)
            else:
                expected_but_not_zero.append(col)
        else:
            if is_all_zero:
                unexpected_zero.append(col)
    
    return unexpected_zero, expected_but_not_zero, expected_and_zero


def _print_results(file_name: str, unexpected_zero: list, expected_but_not_zero: list, expected_and_zero: list) -> bool:
    """Print results and return True if passed."""
    passed = True
    
    if expected_and_zero:
        print(f"INFO: {file_name} verified expected zero columns: {', '.join(expected_and_zero)}")
    
    if expected_but_not_zero:
        print(f"WARNING: {file_name} has expected-zero columns that are NOT zero: {', '.join(expected_but_not_zero)}")
    
    if unexpected_zero:
        print(f"ERROR: {file_name} has unexpected columns with all zero values: {', '.join(unexpected_zero)}")
        passed = False
    
    if passed:
        print(f"PASS: {file_name} has no unexpected columns with all zero values")
    
    return passed


def check_gene_metrics(file_path: str, expected_zero: list = None) -> bool:
    """Check for columns that are all zeros in a gene metrics file."""
    if expected_zero is None:
        expected_zero = GENE_METRICS_EXPECTED_ZERO
    
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return False
    
    try:
        df = _load_csv(file_path)
        unexpected_zero, expected_but_not_zero, expected_and_zero = _find_all_zero_columns(df, expected_zero)
        return _print_results("Gene Metrics", unexpected_zero, expected_but_not_zero, expected_and_zero)
    except Exception as e:
        print(f"ERROR: Could not process Gene Metrics: {str(e)}")
        return False


def check_cell_metrics(file_path: str, expected_zero: list = None) -> bool:
    """Check for columns that are all zeros in a cell metrics file."""
    if expected_zero is None:
        expected_zero = CELL_METRICS_EXPECTED_ZERO
    
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return False
    
    try:
        df = _load_csv(file_path)
        unexpected_zero, expected_but_not_zero, expected_and_zero = _find_all_zero_columns(df, expected_zero)
        return _print_results("Cell Metrics", unexpected_zero, expected_but_not_zero, expected_and_zero)
    except Exception as e:
        print(f"ERROR: Could not process Cell Metrics: {str(e)}")
        return False


def check_library_metrics(file_path: str, expected_zero: list = None) -> bool:
    """Check for metrics with zero values in a library metrics file (key-value format)."""
    if expected_zero is None:
        expected_zero = LIBRARY_METRICS_EXPECTED_ZERO
    
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return False
    
    try:
        # Library metrics are key-value pairs in CSV format
        df = _load_csv(file_path)
        df.columns = ['metric', 'value']
        df['value'] = pd.to_numeric(df['value'], errors='coerce')
        
        unexpected_zero = []
        expected_but_not_zero = []
        expected_and_zero = []
        
        for _, row in df.iterrows():
            if pd.notna(row['value']):
                metric = row['metric']
                is_zero = row['value'] == 0
                
                if metric in expected_zero:
                    if is_zero:
                        expected_and_zero.append(metric)
                    else:
                        expected_but_not_zero.append(metric)
                else:
                    if is_zero:
                        unexpected_zero.append(metric)
        
        return _print_library_results(unexpected_zero, expected_but_not_zero, expected_and_zero)
            
    except Exception as e:
        print(f"ERROR: Could not process Library Metrics: {str(e)}")
        return False


def _print_library_results(unexpected_zero: list, expected_but_not_zero: list, expected_and_zero: list) -> bool:
    """Print results for library metrics and return True if passed."""
    passed = True
    
    if expected_and_zero:
        print(f"INFO: Library Metrics verified expected zero values: {', '.join(expected_and_zero)}")
    
    if expected_but_not_zero:
        print(f"WARNING: Library Metrics has expected-zero metrics that are NOT zero: {', '.join(expected_but_not_zero)}")
    
    if unexpected_zero:
        print(f"ERROR: Library Metrics has unexpected metrics with zero values: {', '.join(unexpected_zero)}")
        passed = False
    
    if passed:
        print(f"PASS: Library Metrics has no unexpected zero values")
    
    return passed


def main():
    parser = argparse.ArgumentParser(
        description='Check for zero metrics in Optimus output files'
    )
    parser.add_argument('--gene-metrics', default=None,
                        help='Path to gene metrics file')
    parser.add_argument('--cell-metrics', default=None,
                        help='Path to cell metrics file')
    parser.add_argument('--library-metrics', default=None,
                        help='Path to library metrics file')
    parser.add_argument('--gene-expected-zero', nargs='*', default=None,
                        help='Gene metrics expected to be zero')
    parser.add_argument('--cell-expected-zero', nargs='*', default=None,
                        help='Cell metrics expected to be zero')
    parser.add_argument('--library-expected-zero', nargs='*', default=None,
                        help='Library metrics expected to be zero')
    
    args = parser.parse_args()
    
    # Check that at least one metrics file is provided
    if not any([args.gene_metrics, args.cell_metrics, args.library_metrics]):
        print("ERROR: At least one metrics file must be provided")
        sys.exit(1)
    
    print("=== This tool currently assumes Optimus Multiome mode, for other modes errors will be shown for expected zero metrics ===")
    print("=== OPTIMUS ZERO METRIC QC RESULTS ===")
    
    results = []
    
    if args.gene_metrics:
        results.append(check_gene_metrics(args.gene_metrics, args.gene_expected_zero))
    else:
        print("SKIPPED: Gene Metrics (not provided)")
    
    if args.cell_metrics:
        results.append(check_cell_metrics(args.cell_metrics, args.cell_expected_zero))
    else:
        print("SKIPPED: Cell Metrics (not provided)")
    
    if args.library_metrics:
        results.append(check_library_metrics(args.library_metrics, args.library_expected_zero))
    else:
        print("SKIPPED: Library Metrics (not provided)")
    
    print("=== OVERALL SUMMARY ===")
    
    if not results:
        print("ERROR: No tests were run")
        sys.exit(1)
    
    if all(results):
        print("OVERALL: All metrics files passed zero value checks")
        sys.exit(0)
    else:
        print("OVERALL: One or more metrics files contain problematic zero values")
        sys.exit(1)


if __name__ == "__main__":
    main()
