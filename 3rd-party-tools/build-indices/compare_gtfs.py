#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import defaultdict

def parse_attributes(attr_str):
    """Parse GTF attributes (9th field) into a dictionary."""
    attrs = {}
    for pair in attr_str.strip().split(';'):
        pair = pair.strip()
        if not pair:
            continue
        key_value = pair.strip().split(' ', 1)
        if len(key_value) == 2:
            key, value = key_value
            attrs[key.strip()] = value.strip('"')
    return attrs

def compare_gtfs(gtf1_path, gtf2_path, output_prefix):
    """Compare two GTF files and generate comparison reports."""
    
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
    
    # Get the minimum length to compare common rows
    min_len = min(len(gtf1), len(gtf2))
    
    for col in cols[:8]:
        # Compare only the common rows
        differences = (gtf1[col].iloc[:min_len] != gtf2[col].iloc[:min_len]).sum()
        if differences > 0:
            struct_diff.append(f"Column '{col}' differs in {differences} rows (of common rows)")
            
            # Sample of differences
            diff_indices = (gtf1[col].iloc[:min_len] != gtf2[col].iloc[:min_len])
            sample_diffs = []
            for idx in diff_indices[diff_indices].index[:5]:  # Show first 5 differences
                sample_diffs.append(f"    Row {idx}: GTF1='{gtf1[col].iloc[idx]}', GTF2='{gtf2[col].iloc[idx]}'")
            struct_diff.extend(sample_diffs)
    
    # Compare attributes (9th field)
    print("\nAnalyzing attribute differences (9th field)...")
    attr_diff = defaultdict(lambda: {'only_in_1': 0, 'only_in_2': 0, 'different_values': 0})
    
    # Only compare the common rows for attributes
    for i in range(min_len):
        attr1 = gtf1['attributes'].iloc[i]
        attr2 = gtf2['attributes'].iloc[i]
        dict1 = parse_attributes(attr1)
        dict2 = parse_attributes(attr2)
        
        # Find all unique keys
        all_keys = set(dict1.keys()) | set(dict2.keys())
        
        for key in all_keys:
            if key in dict1 and key not in dict2:
                attr_diff[key]['only_in_1'] += 1
            elif key in dict2 and key not in dict1:
                attr_diff[key]['only_in_2'] += 1
            elif dict1.get(key) != dict2.get(key):
                attr_diff[key]['different_values'] += 1
    
    # Write reports
    with open(f"{output_prefix}_structural_diff.txt", 'w') as f:
        f.write("Structural Differences (first 8 fields):\n")
        for diff in struct_diff:
            f.write(f"{diff}\n")
    
    with open(f"{output_prefix}_attribute_diff.txt", 'w') as f:
        f.write("Attribute Differences (9th field):\n\n")
        f.write(f"{'Attribute Key':<30} {'Only in GTF1':<15} {'Only in GTF2':<15} {'Different Values':<15}\n")
        f.write("-" * 75 + "\n")
        for key, counts in attr_diff.items():
            if any(counts.values()):  # Only write if there are differences
                f.write(f"{key:<30} {counts['only_in_1']:<15} {counts['only_in_2']:<15} {counts['different_values']:<15}\n")
    
    # Generate gene-level comparison
    print("\nAnalyzing gene-level differences...")
    genes1 = set()
    genes2 = set()
    
    for attr in gtf1['attributes']:
        attrs = parse_attributes(attr)
        if 'gene_id' in attrs:
            genes1.add(attrs['gene_id'])
    
    for attr in gtf2['attributes']:
        attrs = parse_attributes(attr)
        if 'gene_id' in attrs:
            genes2.add(attrs['gene_id'])
    
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
        
        # Add a section for MT genes comparison
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
    
    print(f"Comparing {args.gtf1} and {args.gtf2}")
    compare_gtfs(args.gtf1, args.gtf2, args.output_prefix)
    print(f"\nComparison complete. Check the following files for results:")
    print(f"  {args.output_prefix}_structural_diff.txt")
    print(f"  {args.output_prefix}_attribute_diff.txt")
    print(f"  {args.output_prefix}_gene_diff.txt")

if __name__ == "__main__":
    main()