#!/usr/bin/env python3
"""
Script for checking zero metrics in Optimus output files.

This script checks for columns/metrics that are all zeros in gene metrics,
cell metrics, and library metrics files. Tests are skipped if the corresponding
file argument is not provided.
"""

import argparse
import pandas as pd
import sys
import os


def _parse_numeric(value):
    if value is None:
        return None
    v = str(value).strip()
    if v == "":
        return None
    # allow percentages and thousands separators
    v = v.replace(',', '')
    if v.endswith('%'):
        v = v[:-1]
    try:
        return float(v)
    except ValueError:
        return None


def check_zero_columns(file_path, file_name, skip_metrics=None):
    """Check for columns that are all zeros in a CSV/TSV file"""
    if skip_metrics is None:
        skip_metrics = ["reads_mapped_exonic_as"]
    
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return False
    
    try:
        # Try reading as CSV first, then TSV
        try:
            df = pd.read_csv(file_path)
        except:
            df = pd.read_csv(file_path, sep='\t')
        
        # Find numeric columns that are all zeros
        zero_columns = []
        skipped_columns = []
        for col in df.columns:
            if pd.api.types.is_numeric_dtype(df[col]):
                if (df[col] == 0).all():
                    if col in skip_metrics:
                        skipped_columns.append(col)
                    else:
                        zero_columns.append(col)
        
        if skipped_columns:
            print(f"INFO: {file_name} skipped zero-value columns: {', '.join(skipped_columns)}")
        
        if zero_columns:
            print(f"ERROR: {file_name} has columns with all zero values: {', '.join(zero_columns)}")
            return False
        else:
            print(f"PASS: {file_name} has no problematic columns with all zero values")
            return True
            
    except Exception as e:
        print(f"ERROR: Could not process {file_name}: {str(e)}")
        return False


def check_zero_values_in_library_metrics(file_path, file_name, skip_metrics=None):
    """Check for metrics with zero values in a key-value format file"""
    if skip_metrics is None:
        skip_metrics = [
            "reads_mapped_antisense_to_gene",
            "reads_mapped_confidently_to_intronic_regions",
            "percent_intronic_reads",
        ]
    
    if not os.path.exists(file_path):
        print(f"ERROR: File not found: {file_path}")
        return False
    
    zero_metrics = []
    skipped_metrics = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                # Try simple key/value parsing (colon, tab, comma)
                if ':' in line and ',' not in line:
                    key, value = line.split(':', 1)
                elif '\t' in line:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        key, value = parts[0], parts[1]
                    else:
                        continue
                elif ',' in line:
                    parts = line.split(',')
                    if len(parts) >= 2:
                        key, value = parts[0], parts[1]
                    else:
                        continue
                else:
                    continue
                
                key = key.strip()
                value = value.strip()
                
                numeric_value = _parse_numeric(value)
                if numeric_value is None:
                    continue
                
                if numeric_value == 0.0:
                    if key in skip_metrics:
                        skipped_metrics.append(key)
                    else:
                        zero_metrics.append(key)
        
        if skipped_metrics:
            print(f"INFO: {file_name} skipped zero-value metrics: {', '.join(skipped_metrics)}")
        
        if zero_metrics:
            print(f"ERROR: {file_name} has metrics with zero values: {', '.join(zero_metrics)}")
            return False
        else:
            print(f"PASS: {file_name} has no problematic metrics with zero values")
            return True
            
    except Exception as e:
        print(f"ERROR: Could not process {file_name}: {str(e)}")
        return False


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
    parser.add_argument('--gene-skip-metrics', nargs='*', default=None,
                        help='Gene metrics to skip')
    parser.add_argument('--cell-skip-metrics', nargs='*', default=None,
                        help='Cell metrics to skip')
    parser.add_argument('--library-skip-metrics', nargs='*', default=None,
                        help='Library metrics to skip')
    
    args = parser.parse_args()
    
    # Check that at least one metrics file is provided
    if not any([args.gene_metrics, args.cell_metrics, args.library_metrics]):
        print("ERROR: At least one metrics file must be provided")
        sys.exit(1)
    
    print("=== OPTIMUS ZERO METRIC QC RESULTS ===")
    
    results = []
    
    # Only run tests for provided metrics files
    if args.gene_metrics:
        gene_passed = check_zero_columns(
            args.gene_metrics, "Gene Metrics", args.gene_skip_metrics
        )
        results.append(gene_passed)
    else:
        print("SKIPPED: Gene Metrics (not provided)")
    
    if args.cell_metrics:
        cell_passed = check_zero_columns(
            args.cell_metrics, "Cell Metrics", args.cell_skip_metrics
        )
        results.append(cell_passed)
    else:
        print("SKIPPED: Cell Metrics (not provided)")
    
    if args.library_metrics:
        library_passed = check_zero_values_in_library_metrics(
            args.library_metrics, "Library Metrics", args.library_skip_metrics
        )
        results.append(library_passed)
    else:
        print("SKIPPED: Library Metrics (not provided)")
    
    print("=== OVERALL SUMMARY ===")
    
    if not results:
        print("ERROR: No tests were run")
        sys.exit(1)
    
    all_passed = all(results)
    
    if not all_passed:
        print("OVERALL: One or more metrics files contain problematic zero values")
        sys.exit(1)
    else:
        print("OVERALL: All metrics files passed zero value checks")
        sys.exit(0)


if __name__ == "__main__":
    main()