import argparse
import pandas as pd
import numpy as np
def merge_matrices(summary_file, align_file, cell_reads, counting_mode):
    # Read the whitelist into a set
    print("Reading Aligning features txt file")
    align = pd.read_csv(align_file, sep="\s+", header=None)
    align.columns=("metric", "value", "shard")    
    align_pv=align.pivot(index="shard", columns="metric", values="value")
    print("Reading summary txt file")
    summary = pd.read_csv(summary_file, sep=",", header=None)
    summary.columns=("metric","value","shard")
    summary_pv=summary.pivot(index="shard", columns="metric", values="value")
    merge_pv=pd.merge(summary_pv, align_pv, left_index=True, right_index=True, how='outer')
    merge_pv.reset_index()
    print("Setting n_reads to numeric")
    for x in merge_pv.columns:
        merge_pv[x]=pd.to_numeric(merge_pv[x], errors='coerce')
    
    n_reads=merge_pv["Number of Reads"].sum()
    
    if counting_mode=="sc_rna":
        counting="Gene"
    else:
        counting="GeneFull_Ex50pAS"
    print(counting)
    print("Calcuating metrics")
    merge_pv[f"Reads Mapped to {counting}: Unique {counting}*n_reads"] = merge_pv[f'Reads Mapped to {counting}: Unique {counting}']*merge_pv['Number of Reads']
    Sum_reads_mapped_unique_gene=merge_pv[f"Reads Mapped to {counting}: Unique {counting}*n_reads"].sum()
    Total_reads_mapped_unique_gene=Sum_reads_mapped_unique_gene/n_reads
    merge_pv["Reads Mapped to Genome: Unique*n_reads"] = merge_pv["Reads Mapped to Genome: Unique"]*merge_pv['Number of Reads']
    merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"] = merge_pv["Reads Mapped to Genome: Unique+Multiple"]*merge_pv['Number of Reads']
    merge_pv["Q30 Bases in RNA read*n_reads"] = merge_pv["Q30 Bases in RNA read"]*merge_pv['Number of Reads']
    merge_pv["Q30 Bases in CB+UMI*n_reads"] = merge_pv["Q30 Bases in CB+UMI"]*merge_pv['Number of Reads']
    merge_pv["Reads With Valid Barcodes*n_reads"] = merge_pv["Reads With Valid Barcodes"]*merge_pv['Number of Reads']
    Sequencing_Saturations_Total = 1-(merge_pv["yesUMIs"].sum()/merge_pv["yessubWLmatch_UniqueFeature"].sum())
    Reads_mapped_Genome_unique=merge_pv["Reads Mapped to Genome: Unique*n_reads"].sum()/n_reads
    Reads_mapped_Genome_unique_multi=merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"].sum()/n_reads
    Q30_RNA=merge_pv["Q30 Bases in RNA read*n_reads"].sum()/n_reads
    Q30_CB_UMI=merge_pv["Q30 Bases in CB+UMI*n_reads"].sum()/n_reads
    valid_barcodes=merge_pv["Reads With Valid Barcodes*n_reads"].sum()/n_reads
    
    print("Calculating cell read metrics")
    
    cells=pd.read_csv(cell_reads, sep="\t")
    cells=cells.drop(cells[cells["CB"]=="CBnotInPasslist"].index)
    cells=cells.sort_values(by=['CB','cbMatch'])
    cells=cells.drop_duplicates(subset="CB", keep='last')
    #cells['shard_number']=cells['shard_number'].apply(str)
    #cells["Index"]=cells["CB"]+"-"+cells['shard_number']
    reads_mapped_antisense_to_gene=cells["exonicAS"].sum()+cells["intronicAS"].sum()
    reads_exonic=cells["exonic"].sum()
    reads_mapped_confidently_to_genome=cells["genomeU"].sum()
    reads_mapped_confidently_to_intronic_regions=cells["intronic"].sum()
    reads_mapped_confidently_to_transcriptome=cells["featureU"].sum()
    total_genes_detected=cells["nGenesUnique"].sum()
    data = {"Number_of_reads": [n_reads], 
        "Sequencing_saturation": [Sequencing_Saturations_Total], 
        "Fraction_of_unique_reads_mapped_to_genome": [Reads_mapped_Genome_unique],
        "Fraction_of_unique_and_multiple_reads_mapped_to_genome": [Reads_mapped_Genome_unique_multi],
        "Fraction_of_reads_with_Q30_bases_in_rna": [Q30_RNA],
        "Fraction_of_reads_with_Q30_bases_in_cb_and_umi": [Q30_CB_UMI],
        "Fraction_of_reads_with_valid_barcodes": [valid_barcodes],
        "Reads_mapped_antisense_to_gene": [reads_mapped_antisense_to_gene],
        "Reads_mapped_confidently_exonic": [reads_exonic],
        "Reads_mapped_confidently_to_genome": [reads_mapped_confidently_to_genome],
        "Reads_mapped_confidently_to_intronic_regions": [reads_mapped_confidently_to_intronic_regions],
        "Reads_mapped_confidently_to_transcriptome": [reads_mapped_confidently_to_transcriptome],
        "Total_genes_detected": [total_genes_detected]
        }
    df=pd.DataFrame(data)
    return df

def main():
    parser = argparse.ArgumentParser(description="Count matching DNA barcodes and determine the best matching method.")
    parser.add_argument("summary_file", help="Path to combined summary metrics.")
    parser.add_argument("align_file", help="Path to combined aligner features metrics.")
    parser.add_argument("cell_reads", help="Path to combined cell reads metrics.")
    parser.add_argument("counting_mode", help="Counting mode for STARsolo alignment.")
    parser.add_argument("base_name", help="How to name output files.")

    args = parser.parse_args()

    df= merge_matrices(args.summary_file, args.align_file, args.cell_reads, args.counting_mode)
    df.transpose().to_csv(args.base_name+"_library_metrics.csv", header=None)

if __name__ == "__main__":
    main()        
