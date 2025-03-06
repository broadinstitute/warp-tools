#!/usr/bin/env python3

import argparse
import pandas as pd
import os
from collections import defaultdict

def ensure_output_directory(output_prefix):
    """Ensure the output directory exists."""
    output_dir = os.path.dirname(output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

def verify_input_files(gtf1_path, gtf2_path):
    """Verify input GTF files exist and are readable."""
    if not os.path.exists(gtf1_path):
        raise FileNotFoundError(f"First GTF file not found: {gtf1_path}")
    if not os.path.exists(gtf2_path):
        raise FileNotFoundError(f"Second GTF file not found: {gtf2_path}")

def parse_attributes(attr_str):
    """Parse GTF attributes (9th field) into a dictionary, normalizing the values."""
    attrs = {}
    if not attr_str or not isinstance(attr_str, str):
        return attrs
        
    for pair in attr_str.strip().split(';'):
        pair = pair.strip()
        if not pair:
            continue
        try:
            key_value = pair.strip().split(' ', 1)
            if len(key_value) == 2:
                key, value = key_value
                # Normalize key and value by stripping quotes and whitespace
                key = key.strip()
                value = value.strip().strip('"')
                attrs[key] = value
        except Exception as e:
            print(f"Warning: Could not parse attribute pair: {pair}")
            continue
    return attrs

def get_attribute_differences(dict1, dict2):
    """Get detailed differences between two attribute dictionaries."""
    extra_in_1 = {k: dict1[k] for k in dict1.keys() - dict2.keys()}
    extra_in_2 = {k: dict2[k] for k in dict2.keys() - dict1.keys()}
    different_values = {}
    same_values = {}
    
    # Check common keys for value differences
    common_keys = set(dict1.keys()) & set(dict2.keys())
    for key in common_keys:
        if dict1[key] != dict2[key]:
            different_values[key] = (dict1[key], dict2[key])
        else:
            same_values[key] = dict1[key]
            
    return extra_in_1, extra_in_2, different_values, same_values

def format_attribute_diff(extra_1, extra_2, diff_vals, row_num=None):
    """Format attribute differences for reporting."""
    lines = []
    if row_num is not None:
        lines.append(f"\nRow {row_num}:")
    
    if extra_1:
        lines.append("  Extra attributes in GTF1:")
        for k, v in sorted(extra_1.items()):
            lines.append(f"    {k} \"{v}\"")
    
    if extra_2:
        lines.append("  Extra attributes in GTF2:")
        for k, v in sorted(extra_2.items()):
            lines.append(f"    {k} \"{v}\"")
    
    if diff_vals:
        lines.append("  Different values:")
        for k, (v1, v2) in sorted(diff_vals.items()):
            lines.append(f"    {k}: GTF1=\"{v1}\", GTF2=\"{v2}\"")
    
    return "\n".join(lines)

def compare_gtfs(gtf1_path, gtf2_path, output_prefix):
    """Compare two GTF files and generate comparison reports."""
    
    # Verify inputs and create output directory
    verify_input_files(gtf1_path, gtf2_path)
    ensure_output_directory(output_prefix)
    
    # Read both GTF files
    print("Reading GTF files...")
    cols = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
    gtf1 = pd.read_csv(gtf1_path, sep='\t', comment='#', names=cols)
    gtf2 = pd.read_csv(gtf2_path, sep='\t', comment='#', names=cols)
    
    # Compare first 8 fields
    print("\nAnalyzing structural differences (first 8 fields)...")
    struct_diff = []
    struct_diff.append(f"Total rows in GTF1: {len(gtf1)}")
    struct_diff.append(f"Total rows in GTF2: {len(gtf2)}")
    struct_diff.append(f"Row difference: {abs(len(gtf1) - len(gtf2))}")
    
    min_len = min(len(gtf1), len(gtf2))
    
    for col in cols[:8]:
        try:
            differences = (gtf1[col].iloc[:min_len] != gtf2[col].iloc[:min_len]).sum()
            if differences > 0:
                struct_diff.append(f"Column '{col}' differs in {differences} rows (of common rows)")
                
                diff_indices = (gtf1[col].iloc[:min_len] != gtf2[col].iloc[:min_len])
                sample_diffs = []
                for idx in diff_indices[diff_indices].index[:5]:
                    sample_diffs.append(f"    Row {idx}: GTF1='{gtf1[col].iloc[idx]}', GTF2='{gtf2[col].iloc[idx]}'")
                struct_diff.extend(sample_diffs)
        except Exception as e:
            print(f"Warning: Error comparing column {col}: {str(e)}")
    
    # Compare attributes (9th field)
    print("\nAnalyzing attribute differences (9th field)...")
    attr_diff = defaultdict(lambda: {'only_in_1': 0, 'only_in_2': 0, 'different_values': 0})
    detailed_diffs = []
    
    for i in range(min_len):
        try:
            attr1 = gtf1['attributes'].iloc[i]
            attr2 = gtf2['attributes'].iloc[i]
            
            dict1 = parse_attributes(attr1)
            dict2 = parse_attributes(attr2)
            
            extra_1, extra_2, diff_vals, same_vals = get_attribute_differences(dict1, dict2)
            
            if extra_1 or extra_2 or diff_vals:
                diff_report = format_attribute_diff(extra_1, extra_2, diff_vals, i)
                detailed_diffs.append(diff_report)
            
            for k in extra_1:
                attr_diff[k]['only_in_1'] += 1
            for k in extra_2:
                attr_diff[k]['only_in_2'] += 1
            for k in diff_vals:
                attr_diff[k]['different_values'] += 1
                
        except Exception as e:
            print(f"Warning: Error comparing attributes in row {i}: {str(e)}")
    
    # Generate gene-level comparison
    print("\nAnalyzing gene-level differences...")
    genes1 = set()
    genes2 = set()
    
    for attr in gtf1['attributes']:
        try:
            attrs = parse_attributes(attr)
            if 'gene_id' in attrs:
                genes1.add(attrs['gene_id'])
        except Exception as e:
            print(f"Warning: Error parsing gene_id: {str(e)}")
    
    for attr in gtf2['attributes']:
        try:
            attrs = parse_attributes(attr)
            if 'gene_id' in attrs:
                genes2.add(attrs['gene_id'])
        except Exception as e:
            print(f"Warning: Error parsing gene_id: {str(e)}")
    
    # Write reports
    with open(f"{output_prefix}_structural_diff.txt", 'w') as f:
        f.write("Structural Differences (first 8 fields):\n")
        for diff in struct_diff:
            f.write(f"{diff}\n")
    
    with open(f"{output_prefix}_attribute_diff.txt", 'w') as f:
        f.write("Attribute Differences (9th field):\n\n")
        f.write(f"{'Attribute Key':<30} {'Only in GTF1':<15} {'Only in GTF2':<15} {'Different Values':<15}\n")
        f.write("-" * 75 + "\n")
        for key, counts in sorted(attr_diff.items()):
            if any(counts.values()):
                f.write(f"{key:<30} {counts['only_in_1']:<15} {counts['only_in_2']:<15} {counts['different_values']:<15}\n")
        
        if detailed_diffs:
            f.write("\nDetailed Attribute Differences:\n")
            f.write("-" * 30 + "\n")
            for diff in detailed_diffs[:10]:
                f.write(f"{diff}\n")
            if len(detailed_diffs) > 10:
                f.write(f"\n... and {len(detailed_diffs) - 10} more differences\n")
    
    with open(f"{output_prefix}_gene_diff.txt", 'w') as f:
        f.write("Gene-level Differences:\n\n")
        f.write(f"Total genes in GTF1: {len(genes1)}\n")
        f.write(f"Total genes in GTF2: {len(genes2)}\n")
        f.write(f"Genes only in GTF1: {len(genes1 - genes2)}\n")
        f.write(f"Genes only in GTF2: {len(genes2 - genes1)}\n")
        f.write(f"Genes in both: {len(genes1 & genes2)}\n\n")
        
        f.write("Genes only in GTF1:\n")
        for gene in sorted(genes1 - genes2):
            f.write(f"{gene}\n")
        
        f.write("\nGenes only in GTF2:\n")
        for gene in sorted(genes2 - genes1):
            f.write(f"{gene}\n")
        
        mt_genes1 = {g for g in genes1 if g.startswith("MT-")}
        mt_genes2 = {g for g in genes2 if g.startswith("MT-")}
        
        f.write("\nMitochondrial Genes Comparison:\n")
        f.write(f"MT genes in GTF1: {len(mt_genes1)}\n")
        f.write(f"MT genes in GTF2: {len(mt_genes2)}\n")
        f.write(f"MT genes only in GTF1: {len(mt_genes1 - mt_genes2)}\n")
        f.write(f"MT genes only in GTF2: {len(mt_genes2 - mt_genes1)}\n")
        
        if mt_genes1 - mt_genes2:
            f.write("\nMT genes only in GTF1:\n")
            for gene in sorted(mt_genes1 - mt_genes2):
                f.write(f"{gene}\n")
        
        if mt_genes2 - mt_genes1:
            f.write("\nMT genes only in GTF2:\n")
            for gene in sorted(mt_genes2 - mt_genes1):
                f.write(f"{gene}\n")

def main():
    parser = argparse.ArgumentParser(description="Compare two GTF files and analyze their differences")
    parser.add_argument("gtf1", help="First GTF file (e.g., from old script)")
    parser.add_argument("gtf2", help="Second GTF file (e.g., from new script)")
    parser.add_argument("--output-prefix", "-o", default="gtf_comparison",
                       help="Prefix for output files (default: gtf_comparison)")
    
    args = parser.parse_args()
    
    try:
        print(f"Comparing {args.gtf1} and {args.gtf2}")
        compare_gtfs(args.gtf1, args.gtf2, args.output_prefix)
        print(f"\nComparison complete. Check the following files for results:")
        print(f"  {args.output_prefix}_structural_diff.txt")
        print(f"  {args.output_prefix}_attribute_diff.txt")
        print(f"  {args.output_prefix}_gene_diff.txt")
    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == "__main__":
    main()