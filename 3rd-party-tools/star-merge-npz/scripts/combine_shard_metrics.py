import argparse
import pandas as pd
import numpy as np

def merge_matrices(summary_file, align_file, cell_reads, counting_mode, uniform_barcodes, uniform_mtx, expected_cells):
    # Read the whitelist into a set
    expected_cells = int(expected_cells)
    print("Reading Aligning features txt file")
    align = pd.read_csv(align_file, sep="\s+", header=None)
    align.columns = ("metric", "value", "shard")
    align_pv = align.pivot(index="shard", columns="metric", values="value")
    
    print("Reading summary txt file")
    summary = pd.read_csv(summary_file, sep=",", header=None)
    summary.columns = ("metric", "value", "shard")
    summary_pv = summary.pivot(index="shard", columns="metric", values="value")
    
    merge_pv = pd.merge(summary_pv, align_pv, left_index=True, right_index=True, how='outer')
    merge_pv.reset_index()
    
    print("Setting n_reads to numeric")
    for x in merge_pv.columns:
        merge_pv[x] = pd.to_numeric(merge_pv[x], errors='coerce')
    
    n_reads = merge_pv["Number of Reads"].sum()
    
    if counting_mode == "sc_rna":
        counting = "Gene"
    else:
        counting = "GeneFull_Ex50pAS"
    
    print(counting)
    
    print("Calculating metrics")
    merge_pv[f"Reads Mapped to {counting}: Unique {counting}*n_reads"] = merge_pv[f'Reads Mapped to {counting}: Unique {counting}'] * merge_pv['Number of Reads']
    sum_reads_mapped_unique_gene = merge_pv[f"Reads Mapped to {counting}: Unique {counting}*n_reads"].sum()
    total_reads_mapped_unique_gene = sum_reads_mapped_unique_gene / n_reads 
    merge_pv["Reads Mapped to Genome: Unique*n_reads"] = merge_pv["Reads Mapped to Genome: Unique"] * merge_pv['Number of Reads']
    merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"] = merge_pv["Reads Mapped to Genome: Unique+Multiple"] * merge_pv['Number of Reads']
    merge_pv["Q30 Bases in RNA read*n_reads"] = merge_pv["Q30 Bases in RNA read"] * merge_pv['Number of Reads']
    merge_pv["Q30 Bases in CB+UMI*n_reads"] = merge_pv["Q30 Bases in CB+UMI"] * merge_pv['Number of Reads']
    merge_pv["Reads With Valid Barcodes*n_reads"] = merge_pv["Reads With Valid Barcodes"] * merge_pv['Number of Reads']  
    sequencing_saturations_total = 1 - (merge_pv["yesUMIs"].sum() / merge_pv["yessubWLmatch_UniqueFeature"].sum())
    reads_mapped_genome_unique = merge_pv["Reads Mapped to Genome: Unique*n_reads"].sum() / n_reads
    reads_mapped_genome_unique_multi = merge_pv["Reads Mapped to Genome: Unique+Multiple*n_reads"].sum() / n_reads
    q30_rna = merge_pv["Q30 Bases in RNA read*n_reads"].sum() / n_reads
    q30_cb_umi = merge_pv["Q30 Bases in CB+UMI*n_reads"].sum() / n_reads
    valid_barcodes = merge_pv["Reads With Valid Barcodes*n_reads"].sum() / n_reads
    
    print("Calculating cell read metrics")
    
    cells = pd.read_csv(cell_reads, sep="\t")
    cells = cells.drop(cells[cells["CB"] == "CBnotInPasslist"].index)
    cells = cells.sort_values(by=['CB', 'cbMatch'])
    cells = cells.drop_duplicates(subset="CB", keep='last')
    
    reads_mapped_antisense_to_gene = cells["exonicAS"].sum() + cells["intronicAS"].sum()
    reads_exonic = cells["exonic"].sum()
    reads_mapped_confidently_to_genome = cells["genomeU"].sum()
    reads_mapped_confidently_to_intronic_regions = cells["intronic"].sum()
    reads_mapped_confidently_to_transcriptome = cells["featureU"].sum()

    print("Reading in combined, filtered matrix market file")
    filtered = pd.read_csv(uniform_mtx, sep="\s+", skiprows=3, header=None)
    unique_columns = filtered[1].unique()
    estimated_cells = len(unique_columns)
    print("Estimated cells are: ", estimated_cells)
    
    umis_in_cells = filtered[2].sum()
    mean_umi_per_cell = umis_in_cells / estimated_cells
    sum_counts_barcode = filtered.groupby(1)[2].sum().reset_index()
    median_umi_per_cell = sum_counts_barcode[2].median()

    print("Reading in uniform filtered barcodes TSV")
    barcodes = pd.read_csv(uniform_barcodes, sep="\t", header=None)
    barcodes.set_index(0, inplace=True)

    print("Subsetting cell_reads to filtered barcodes")
    cells.set_index("CB", inplace=True)
    cells_filtered = cells[cells.index.isin(barcodes.index)]

    print("Length of cells_filtered: ", len(cells_filtered))
    print("Estimated cells (should match above): ", estimated_cells)

    unique_reads_in_cells_mapped_to_gene = cells_filtered["countedU"].sum()
    fraction_of_unique_reads_in_cells = unique_reads_in_cells_mapped_to_gene / sum_reads_mapped_unique_gene
    mean_reads_per_cell = unique_reads_in_cells_mapped_to_gene/len(cells_filtered)
    median_reads_per_cell = cells_filtered["countedU"].median()
    mean_gene_per_cell = cells_filtered["nGenesUnique"].sum()/len(cells_filtered)
    median_gene_per_cell = cells_filtered["nGenesUnique"].median()
    unique_rows = filtered[0].unique()
    total_genes_unique_detected = len(unique_rows)

    # Keeper cell metrics
    # Expected number of cells is unknown, so this commented out for now
    #expected_cells = 3000  # Placeholder, replace with actual value
    percent_target = estimated_cells/expected_cells
    percent_intronic_reads = reads_mapped_confidently_to_intronic_regions/n_reads

    if counting_mode == "sc_rna":
        gene_threshold = 1500
    else:
        gene_threshold = 1000

    keeper_cells = cells_filtered[cells_filtered["nGenesUnique"] > gene_threshold]
    keeper_mean_reads_per_cell = keeper_cells["countedU"].mean()
    keeper_median_genes = keeper_cells["nGenesUnique"].median()
    keeper_cells_count = len(keeper_cells)
    percent_keeper = keeper_cells_count/estimated_cells
    percent_usable = keeper_cells_count/expected_cells

    data = {
        "number_of_reads": [n_reads],
        "sequencing_saturation": [sequencing_saturations_total],
        "fraction_of_unique_reads_mapped_to_genome": [reads_mapped_genome_unique],
        "fraction_of_unique_and_multiple_reads_mapped_to_genome": [reads_mapped_genome_unique_multi],
        "fraction_of_reads_with_Q30_bases_in_rna": [q30_rna],
        "fraction_of_reads_with_Q30_bases_in_cb_and_umi": [q30_cb_umi],
        "fraction_of_reads_with_valid_barcodes": [valid_barcodes],
        "reads_mapped_antisense_to_gene": [reads_mapped_antisense_to_gene],
        "reads_mapped_confidently_exonic": [reads_exonic],
        "reads_mapped_confidently_to_genome": [reads_mapped_confidently_to_genome],
        "reads_mapped_confidently_to_intronic_regions": [reads_mapped_confidently_to_intronic_regions],
        "reads_mapped_confidently_to_transcriptome": [reads_mapped_confidently_to_transcriptome],
        "estimated_cells": [estimated_cells],
        "umis_in_cells": [umis_in_cells],
        "mean_umi_per_cell": [mean_umi_per_cell],
        "median_umi_per_cell": [median_umi_per_cell],
        "unique_reads_in_cells_mapped_to_gene": [unique_reads_in_cells_mapped_to_gene],
        "fraction_of_unique_reads_in_cells": [fraction_of_unique_reads_in_cells],
        "mean_reads_per_cell": [mean_reads_per_cell],
        "median_reads_per_cell": [median_reads_per_cell],
        "mean_gene_per_cell": [mean_gene_per_cell],
        "median_gene_per_cell": [median_gene_per_cell],
        "total_genes_unique_detected": [total_genes_unique_detected],
        "percent_target": [percent_target],
        "percent_intronic_reads": [percent_intronic_reads],
        "keeper_mean_reads_per_cell": [keeper_mean_reads_per_cell],
        "keeper_median_genes": [keeper_median_genes],
        "keeper_cells": [keeper_cells_count],
        "percent_keeper": [percent_keeper],
        "percent_usable": [percent_usable]
    }

    df = pd.DataFrame(data)
    return df

def main():
    parser = argparse.ArgumentParser(description="Count matching DNA barcodes and determine the best matching method.")
    parser.add_argument("summary_file", help="Path to combined summary metrics.")
    parser.add_argument("align_file", help="Path to combined aligner features metrics.")
    parser.add_argument("cell_reads", help="Path to combined cell reads metrics.")
    parser.add_argument("counting_mode", help="Counting mode for STARsolo alignment.")
    parser.add_argument("base_name", help="How to name output files.")
    parser.add_argument("uniform_barcodes", help="Barcodes TSV returned from STARsolo filtered matrix")
    parser.add_argument("uniform_mtx", help="Matrix Mtx file returned from STARsolo filtered matrix")
    parser.add_argument("expected_cells", help="Expected number of cells based on the number loaded in experiment; default is usually 3000")

    args = parser.parse_args()

    df= merge_matrices(args.summary_file, args.align_file, args.cell_reads, args.counting_mode, args.uniform_barcodes, args.uniform_mtx, args.expected_cells)
    df.transpose().to_csv(args.base_name+"_library_metrics.csv", header=None)

if __name__ == "__main__":
    main()    