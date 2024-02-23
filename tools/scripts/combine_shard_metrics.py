import argparse
import pandas as pd
import numpy as np
def merge_matrices(summary_file, align_file, counting_mode, base_name):
    # Read the whitelist into a set
    align = pd.read_csv("example_align_features.txt", sep="\s+", header=None)
    align.columns=("metric", "value", "shard")    
    align_pv=align.pivot(index="shard", columns="metric", values="value")
    summary = pd.read_csv("example_summary.txt", sep=",", header=None)
    summary.columns=("metric","value","shard")
    summary_pv=summary.pivot(index="shard", columns="metric", values="value")
    merge_pv=pd.merge(summary_pv, align_pv, left_index=True, right_index=True, how='outer')
    merge_pv.reset_index()
    for x in merge_pv.columns:
        merge_pv[x]=pd.to_numeric(merge_pv[x], errors='coerce')
    n_reads=merge_pv["Number of Reads"].sum()
    
    if counting_mode=="sc_rna":
        counting="Gene"
    else:
      counting="GeneFull_Ex50pAS"
    merge_pv[f"Reads Mapped to {counting}: Unique {counting}"]=pd.to_numeric(merge_pv[f"Reads Mapped to {counting}: Unique {counting}"])
    Sum_reads_mapped_unique_gene=merge_pv[f"Reads Mapped to {counting}: Unique {counting}*n_reads"].sum()
    Total_reads_mapped_unique_gene=Sum_reads_mapped_unique_gene/n_reads
    merge_pv["Reads Mapped to Genome: Unique*n_reads"] = merge_pv["Reads Mapped to Genome: Unique"]*merge_pv['Number of Reads']
    merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"] = merge_pv["Reads Mapped to Genome: Unique+Multiple"]*merge_pv['Number of Reads']
    merge_pv["Q30 Bases in RNA read*n_reads"] = merge_pv["Q30 Bases in RNA read"]*merge_pv['Number of Reads']
    merge_pv["Q30 Bases in CB+UMI*n_reads"] = merge_pv["Q30 Bases in CB+UMI"]*merge_pv['Number of Reads']
    merge_pv["Reads With Valid Barcodes*n_reads"] = merge_pv["Reads With Valid Barcodes"]*merge_pv['Number of Reads']
    Sequencing_Saturations_Total = 1-(merge_pv["yesUMIs"].sum()/merge_pv["yessubWLmatch_UniqueFeature"].sum())
    Reads_mapped_Genome_unique=merge_pv["Reads Mapped to Genome: Unique*n_reads"].sum()/n_reads
    Reads_mapped_Genome_unique_multu=merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"].sum()/n_reads
    Q30_RNA=merge_pv["Q30 Bases in RNA read*n_reads"].sum()/n_reads
    Q30_CB_UMI=merge_pv["Q30 Bases in CB+UMI*n_reads"].sum()/n_reads
    valid_barcodes=merge_pv["Reads With Valid Barcodes*n_reads"].sum()/n_reads
    data = {"Number of Reads": [n_reads], 
        "Sequencing Saturation": [Sequencing_Saturations_Total], 
        "Fraction of Unique Reads Mapped to Genome": [Reads_mapped_Genome_unique],
        "Fraction of Unique and Multiple Reads Mapped to Genome": [Reads_mapped_Genome_unique_multu],
        "Fraction of reads with Q30 Bases in RNA": [Q30_RNA],
        "Fraction of reads with Q30 Bases in CB and UMI": [Q30_CB_UMI],
        "Fraction of Reads with Valid Barcodes": [valid_barcodes]}
    df=pd.DataFrame(data)
    return df

df= merge_matrices()
df.transpose().to_csv("test.csv", header=None)

def main():
    parser = argparse.ArgumentParser(description="Count matching DNA barcodes and determine the best matching method.")
    parser.add_argument("summary_file", help="Path to the input file containing DNA barcodes.")
    parser.add_argument("align_file", help="Path to the whitelist file containing allowed barcodes.")
    parser.add_argument("counting_mode", help="Path to the output file for the best matching method parameter.")
    parser.add_argument("base_name", help="Path to the output file for the best matching method parameter.")


    args = parser.parse_args()
    # Define input and output files here you remove argparse
    #input_file = "downsample.fq"  # Change to the path of your input file
    #whitelist_file = "~{whitelist}"  # Change to the path of your whitelist file
    #output_file = "best_match.txt"  # Change to the path of your output file

    (
        count_forward_start,
        count_reverse_comp_start,
        count_forward_end,
        count_reverse_comp_end
    ) = count_matching_barcodes(args.input_file, args.whitelist_file)

if __name__ == "__main__":
    main()        
