#!/usr/bin/env python3
"""
mito_plot_filter.py

Runs mitochondrial post-processing on a Hail MatrixTable:
  - Exports a filtered VCF (bgzipped + tabix-indexed)
  - Exports a sample metadata TSV
  - Generates 7 QC/summary SVG plots
  - Copies all outputs to a cloud destination

Usage:
    python3 mito_plot_filter.py \
        --input-path  gs://bucket/path/to/mt \
        --output-path gs://bucket/output/dir \
        --output-base my_run_prefix
"""

import argparse
import os
import hail as hl
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def parse_args():
    parser = argparse.ArgumentParser(description="Mito post-processing: filter, annotate, and plot.")
    parser.add_argument("--input-path",  required=True, help="GCS path to input Hail MatrixTable")
    parser.add_argument("--output-path", required=True, help="GCS output directory (no trailing slash required)")
    parser.add_argument("--output-base", required=True, help="Prefix for all output file names")
    return parser.parse_args()


def main():
    args = parse_args()

    input_matrix_path = args.input_path
    output_path       = args.output_path
    output_base       = args.output_base

    cloud_prefix = f"{output_path.rstrip('/')}/{output_base}"

    print(f"[mito_plot_filter] output_base  : {output_base}")
    print(f"[mito_plot_filter] cloud_prefix : {cloud_prefix}")
    print(f"[mito_plot_filter] cwd          : {os.getcwd()}")

    hl.init(default_reference="GRCh38", idempotent=True)

    filt_annotated_mt_500k = hl.read_matrix_table(input_matrix_path)

    # -------------------------------------------------------------------------
    # Local output file names (will be copied to cloud at the end)
    # -------------------------------------------------------------------------
    vcf_local             = f"{output_base}.vcf.bgz"
    sample_metadata_local = f"{output_base}_metadata.tsv"

    variants_per_sample_svg            = f"{output_base}.variants_per_sample.svg"
    mito_cn_distribution_svg           = f"{output_base}.mito_cn_distribution.svg"
    variant_allele_frequency_svg       = f"{output_base}.variant_allele_frequency.svg"
    variant_af_and_allele_fraction_svg = f"{output_base}.variant_af_and_allele_fraction.svg"
    numt_fp_by_mtcn_svg                = f"{output_base}.numt_fp_by_mtcn.svg"
    haplogroup_heteroplasmy_svg        = f"{output_base}.haplogroup_heteroplasmy.svg"
    haplogroup_homoplasmy_svg          = f"{output_base}.haplogroup_homoplasmy.svg"

    mt = filt_annotated_mt_500k

    # -------------------------------------------------------------------------
    # Plot 1: Histogram of variants per sample
    # -------------------------------------------------------------------------
    print("Making plot 1")
    mt_tmp = mt.annotate_cols(
        n_variants=hl.agg.count_where(mt.AD[1] > 0)
    )
    df = mt_tmp.cols().select("n_variants").to_pandas()
    plt.figure(figsize=(10, 6))
    sns.histplot(data=df, x="n_variants", bins=50, kde=False, color="#4C72B0")
    plt.title("Distribution of Variants per Sample (Based on AD > 0)", fontsize=15)
    plt.xlabel("Number of Variants (Alt Depth > 0)", fontsize=12)
    plt.ylabel("Count of Samples", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.tight_layout()
    plt.savefig(variants_per_sample_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # Plot 2: Mitochondrial copy number distribution
    # -------------------------------------------------------------------------
    df = mt.cols().select("mito_cn").to_pandas().dropna()
    plt.figure(figsize=(10, 6))
    sns.histplot(df["mito_cn"], bins=50, color="#55A868")
    plt.title("Distribution of Mitochondrial Copy Number", fontsize=15)
    plt.xlabel("Mitochondrial Copy Number", fontsize=12)
    plt.ylabel("Count of Samples", fontsize=12)
    median_cn = df["mito_cn"].median()
    plt.axvline(median_cn, color="red", linestyle="--", label=f"Median: {median_cn:.1f}")
    plt.legend()
    plt.tight_layout()
    plt.savefig(mito_cn_distribution_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # Plot 3: Allele Frequency using AC/AN AF
    # -------------------------------------------------------------------------
    mt_tmp = mt.annotate_rows(
        calculated_AF=(mt.AC_het + mt.AC_hom) / mt.AN
    )
    mt_tmp = mt_tmp.filter_rows(mt_tmp.AN > 0)
    df = mt_tmp.rows().select("calculated_AF").to_pandas()
    plt.figure(figsize=(10, 6))
    sns.histplot(df["calculated_AF"], bins=50, color="#8172B3", log_scale=(False, True))
    plt.title("Distribution of Variant Allele Frequency (Population Level)", fontsize=15)
    plt.xlabel("Allele Frequency (AF)", fontsize=12)
    plt.ylabel("Count of Variants (Log Scale)", fontsize=12)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(variant_allele_frequency_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # Plot 4: AF and allele fraction (side-by-side histograms)
    # -------------------------------------------------------------------------
    mt_tmp = mt.annotate_rows(
        calculated_AF=(mt.AC_het + mt.AC_hom) / mt.AN
    )
    mt_tmp = mt_tmp.filter_rows(mt_tmp.AN > 0)
    af_hist = mt_tmp.aggregate_rows(hl.agg.hist(mt_tmp.calculated_AF, 0, 1, 50))
    mt_tmp = mt_tmp.annotate_entries(
        allele_fraction=hl.cond(
            (mt_tmp.AD[0] + mt_tmp.AD[1]) > 0,
            mt_tmp.AD[1] / (mt_tmp.AD[0] + mt_tmp.AD[1]),
            hl.null(hl.tfloat)
        )
    )
    frac_hist = mt_tmp.aggregate_entries(hl.agg.hist(mt_tmp.allele_fraction, 0, 1, 50))

    fig, ax = plt.subplots(1, 2, figsize=(16, 6))
    af_bin_edges = af_hist.bin_edges
    af_bin_freq  = af_hist.bin_freq
    af_bin_centers = [(af_bin_edges[i] + af_bin_edges[i + 1]) / 2 for i in range(len(af_bin_edges) - 1)]
    ax[0].bar(
        af_bin_centers,
        af_bin_freq,
        width=(af_bin_edges[1] - af_bin_edges[0]),
        color="#8172B3",
        alpha=0.7,
        edgecolor="black"
    )
    ax[0].set_yscale("log")
    ax[0].set_title("Variant Allele Frequency (Population-Level)", fontsize=15)
    ax[0].set_xlabel("Allele Frequency")
    ax[0].set_ylabel("Count (log)")
    ax[0].grid(axis="y", linestyle="--", alpha=0.5)

    frac_bin_edges   = frac_hist.bin_edges
    frac_bin_freq    = frac_hist.bin_freq
    frac_bin_centers = [(frac_bin_edges[i] + frac_bin_edges[i + 1]) / 2 for i in range(len(frac_bin_edges) - 1)]
    ax[1].bar(
        frac_bin_centers,
        frac_bin_freq,
        width=(frac_bin_edges[1] - frac_bin_edges[0]),
        color="#55A868",
        alpha=0.7,
        edgecolor="black"
    )
    ax[1].set_title("Allele Fraction (Per Sample Heteroplasmy)", fontsize=15)
    ax[1].set_xlabel("Allele Fraction (AD1 / DP)")
    ax[1].set_ylabel("Count")
    ax[1].set_xlim(0, 1)
    ax[1].grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(variant_af_and_allele_fraction_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # NUMT flagging + mtCN bins  (used for Plot 5 and later VCF export)
    # -------------------------------------------------------------------------
    mt = filt_annotated_mt_500k
    mt = mt.annotate_entries(
        allele_fraction=hl.cond(
            hl.sum(mt.AD) > 0,
            mt.AD[1] / hl.sum(mt.AD),
            hl.null(hl.tfloat)
        )
    )
    het_expr  = (
        (mt.allele_fraction >= 0.01) &
        (mt.allele_fraction <= 0.05)
    )
    numt_expr = het_expr & mt.common_low_heteroplasmy
    mt = mt.annotate_cols(
        num_hets   =hl.agg.count_where(het_expr & (~mt.common_low_heteroplasmy)),
        num_numt_fp=hl.agg.count_where(numt_expr)
    )
    mt = mt.annotate_cols(
        mtcn_bin=(
            hl.case()
            .when(mt.mito_cn < 25,  "0–25")
            .when(mt.mito_cn < 50,  "25–50")
            .when(mt.mito_cn < 75,  "50–75")
            .when(mt.mito_cn < 100, "75–100")
            .when(mt.mito_cn < 125, "100–125")
            .when(mt.mito_cn < 150, "125–150")
            .when(mt.mito_cn < 175, "150–175")
            .when(mt.mito_cn < 200, "175–200")
            .when(mt.mito_cn < 225, "200–225")
            .when(mt.mito_cn < 250, "225–250")
            .when(mt.mito_cn < 275, "250–275")
            .when(mt.mito_cn < 300, "275–300")
            .default(">=300")
        )
    )
    mt = mt.annotate_cols(
        numt_fp_risk_tier=(
            hl.case()
            .when(mt.mito_cn < 50,  "high")
            .when(mt.mito_cn < 100, "elevated")
            .default("standard")
        )
    )

    ht_samples = mt.cols().select("mtcn_bin", "num_hets", "num_numt_fp")
    df = ht_samples.to_pandas()
    df = df.dropna(subset=["mtcn_bin"])
    bin_order = [
        "50–75", "75–100", "100–125", "125–150", "150–175",
        "175–200", "200–225", "225–250", "250–275", "275–300", ">=300"
    ]
    df_B = (
        df.groupby("mtcn_bin")["num_numt_fp"]
        .mean()
        .reset_index(name="mean_numt_fp")
    )
    df_B["mtcn_bin"] = pd.Categorical(df_B["mtcn_bin"], categories=bin_order, ordered=True)
    df_B = df_B.sort_values("mtcn_bin")

    # -------------------------------------------------------------------------
    # Plot 5: NUMT false positives by mtCN
    # -------------------------------------------------------------------------
    plt.figure(figsize=(8, 5))
    sns.lineplot(data=df_B, x="mtcn_bin", y="mean_numt_fp", marker="o")
    plt.title("Mean NUMT False Positives (0.01–0.50) by mtCN")
    plt.xlabel("mtDNA Copy Number Bin")
    plt.ylabel("Mean Number of NUMT-FPs")
    plt.tight_layout()
    plt.savefig(numt_fp_by_mtcn_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # Haplogroup plots
    # -------------------------------------------------------------------------
    print("Making haplogroup plots")
    mt = filt_annotated_mt_500k
    mt = mt.annotate_entries(
        VAF=hl.if_else(
            hl.len(mt.AD) > 1,
            mt.AD[1] / hl.sum(mt.AD),
            hl.null(hl.tfloat)
        )
    )
    mt = mt.annotate_cols(
        macro_haplogroup=mt.major_haplogroup[0]
    )
    mt = mt.annotate_entries(
        is_heteroplasmic=((mt.VAF >= 0.10) & (mt.VAF <= 0.95))
    )
    mt = mt.annotate_cols(
        n_het_snvs=hl.agg.count_where(mt.is_heteroplasmic)
    )
    df = (
        mt.cols()
        .select("hap", "n_het_snvs")
        .to_pandas()
        .dropna(subset=["hap"])
    )

    def haplogroup_origin(hg):
        if hg.startswith("L"):
            return "African"
        elif hg.startswith(("M", "D", "C", "A", "B", "F")):
            return "Asian"
        return "European"

    df["origin"] = df["hap"].apply(haplogroup_origin)
    origin_colors = {"African": "purple", "Asian": "green", "European": "blue"}
    hap_order = sorted(df["hap"].unique())

    # Plot 6: Heteroplasmic variants by haplogroup
    plt.figure(figsize=(14, 5))
    ax = sns.boxplot(
        data=df,
        x="hap",
        y="n_het_snvs",
        hue="origin",
        order=hap_order,
        palette=origin_colors,
        dodge=False,
        showfliers=False,
        linewidth=0.8
    )
    ax.set_xlabel("mtDNA haplogroup")
    ax.set_ylabel("Heteroplasmic SNVs per sample (VAF 0.10–0.95)")
    ax.set_title("Heteroplasmic variants by mtDNA haplogroup")
    plt.xticks(rotation=90)
    plt.legend(title="Phylogenetic origin")
    plt.tight_layout()
    plt.savefig(haplogroup_heteroplasmy_svg, format="svg")
    plt.close()

    mt = mt.annotate_entries(
        is_homoplasmic=mt.VAF >= 0.95
    )
    mt = mt.annotate_cols(
        n_hom_snvs=hl.agg.count_where(mt.is_homoplasmic)
    )
    df_hom = (
        mt.cols()
        .select("hap", "n_hom_snvs")
        .to_pandas()
        .dropna(subset=["hap"])
    )
    df_hom["origin"] = df_hom["hap"].apply(haplogroup_origin)
    hap_order = sorted(df["hap"].unique())

    # Plot 7: Homoplasmic variants by haplogroup
    plt.figure(figsize=(8, 5))
    ax = sns.boxplot(
        data=df_hom,
        x="hap",
        y="n_hom_snvs",
        order=hap_order,
        hue="origin",
        showfliers=False,
        linewidth=0.8
    )
    ax.set_xlabel("mtDNA hap")
    ax.set_ylabel("Homoplasmic SNVs per sample (VAF ≥ 0.95)")
    ax.set_title("Homoplasmic variant burden by mtDNA haplogroup")
    plt.tight_layout()
    plt.savefig(haplogroup_homoplasmy_svg, format="svg")
    plt.close()

    # -------------------------------------------------------------------------
    # Add NUMT warning entry field and export filtered VCF
    # -------------------------------------------------------------------------
    mt = filt_annotated_mt_500k
    mt = mt.annotate_entries(
        allele_fraction=hl.cond(
            hl.sum(mt.AD) > 0,
            mt.AD[1] / hl.sum(mt.AD),
            hl.null(hl.tfloat)
        )
    )
    het_expr  = (
        (mt.allele_fraction >= 0.01) &
        (mt.allele_fraction <= 0.05)
    )
    numt_expr = het_expr & mt.common_low_heteroplasmy
    mt = mt.annotate_cols(
        num_hets   =hl.agg.count_where(het_expr & (~mt.common_low_heteroplasmy)),
        num_numt_fp=hl.agg.count_where(numt_expr)
    )
    mt = mt.annotate_cols(
        numt_fp_risk_tier=(
            hl.case()
            .when(mt.mito_cn < 50,  "high")
            .when(mt.mito_cn < 100, "elevated")
            .default("standard")
        ),
        low_mtcn_numt_risk=mt.mito_cn < 100
    )
    mt = mt.annotate_entries(
        numt_fp_warning=(
            (mt.allele_fraction >= 0.01) &
            (mt.allele_fraction < 0.05) &
            mt.common_low_heteroplasmy &
            mt.low_mtcn_numt_risk
        )
    )

    mt_vcf = mt.annotate_rows(
        rsid=hl.if_else(
            hl.is_defined(mt.rsid) & (hl.len(mt.rsid) > 0),
            hl.delimit(hl.array(mt.rsid), ";"),
            hl.missing(hl.tstr)
        )
    )
    mt_vcf = mt_vcf.filter_rows(~mt_vcf.filters.contains("npg"))
    
    mt_vcf = mt_vcf.select_rows(
        rsid=mt_vcf.rsid,
        filters=mt_vcf.filters,
        AN=mt_vcf.AN,
        AC_hom=mt_vcf.AC_hom,
        AC_het=mt_vcf.AC_het,
        AF_hom=mt_vcf.AF_hom,
        AF_het=mt_vcf.AF_het,
        dp_mean=mt_vcf.dp_mean,
        mq_mean=mt_vcf.mq_mean,
        tlod_mean=mt_vcf.tlod_mean,
        max_hl=mt_vcf.max_hl,
        hap_defining_variant=mt_vcf.hap_defining_variant,
        variant_context=mt_vcf.variant_context,
        common_low_heteroplasmy=mt_vcf.common_low_heteroplasmy
    )
    mt_vcf = mt_vcf.select_entries(
        GT=mt_vcf.GT,
        DP=mt_vcf.DP,
        AD=mt_vcf.AD,
        HL=mt_vcf.HL,
        TLOD=mt_vcf.TLOD,
        FT=mt_vcf.FT,
        numt_fp_warning=hl.int32(mt_vcf.numt_fp_warning)
    )
    hl.export_vcf(mt_vcf, vcf_local, tabix=True)

    # -------------------------------------------------------------------------
    # Export sample metadata TSV
    # -------------------------------------------------------------------------
    cols_ht = mt.cols()
    sample_ht = cols_ht.select(
        participant_id     =cols_ht.participant_id,
        mito_cn            =cols_ht.mito_cn,
        mt_mean_coverage   =cols_ht.mt_mean_coverage,
        wgs_median_coverage=cols_ht.wgs_median_coverage,
        major_haplogroup   =cols_ht.major_haplogroup,
        contamination      =cols_ht.contamination,
        freemix_percentage =cols_ht.freemix_percentage,
        contam_high_het    =cols_ht.contam_high_het,
        numt_fp_risk_tier  =cols_ht.numt_fp_risk_tier
    )
    sample_ht.export(sample_metadata_local)


if __name__ == "__main__":
    main()
