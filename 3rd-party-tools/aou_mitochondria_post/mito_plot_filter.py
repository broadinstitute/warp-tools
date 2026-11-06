#!/usr/bin/env python3
"""
mito_plot_filter.py

Runs mitochondrial post-processing on a Hail MatrixTable:
    - Adds NUMT warning annotations to entries
  - Exports a filtered VCF (bgzipped + tabix-indexed)
  - Exports a sample metadata TSV

Usage:
    python3 mito_plot_filter.py \
        --input-path  gs://bucket/path/to/mt \
        --output-root gs://bucket/output/dir \
        --basename my_run_prefix
"""

import argparse
import os
import hail as hl


def parse_args():
    parser = argparse.ArgumentParser(description="Mito post-processing: add NUMT warning and export outputs.")
    parser.add_argument("--input-path",  required=True, help="GCS path to input Hail MatrixTable")
    parser.add_argument("--output-root", required=True, help="Cloud output directory root (no trailing slash required)")
    parser.add_argument("--basename", required=True, help="Shared base name used to build all output file names")
    parser.add_argument("--tmp-dir", default="./hail_tmp",
                        help="Local directory for Hail temporary files (should be on the provisioned disk, not /tmp)")
    return parser.parse_args()


def main():
    args = parse_args()

    input_matrix_path = args.input_path
    output_root       = args.output_root
    basename          = args.basename
    tmp_dir           = args.tmp_dir
    output_prefix     = f"{output_root.rstrip('/')}/{basename}"

    print(f"[mito_plot_filter] basename     : {basename}")
    print(f"[mito_plot_filter] output_root  : {output_root}")
    print(f"[mito_plot_filter] output_prefix: {output_prefix}")
    print(f"[mito_plot_filter] tmp_dir      : {tmp_dir}")
    print(f"[mito_plot_filter] cwd          : {os.getcwd()}")

    os.makedirs(tmp_dir, exist_ok=True)

    # Initialize Hail once in idempotent mode so reruns in the same process are safe.
    # tmp_dir is set to a path on the provisioned execution disk (not /tmp) to avoid
    # filling the OS boot disk with Hail spill files when running in Cromwell/Terra WDL.
    hl.init(default_reference="GRCh38", idempotent=True, tmp_dir=tmp_dir)

    # Input is the filtered/annotated MatrixTable produced earlier in the WDL.
    mt = hl.read_matrix_table(input_matrix_path)

    # -------------------------------------------------------------------------
    # Output file paths; these can be local paths or cloud URIs (for example gs://...).
    # -------------------------------------------------------------------------
    vcf_output_path             = f"{output_prefix}.filtered.vcf.bgz"
    sample_metadata_output_path = f"{output_prefix}.metadata.tsv"


    # -------------------------------------------------------------------------
    # Add NUMT warning entry field (matches notebook logic)
    # -------------------------------------------------------------------------
    # Recompute per-entry allele fraction from AD on the split MT.
    mt = mt.annotate_entries(
        allele_fraction=hl.cond(
            hl.sum(mt.AD) > 0,
            mt.AD[1] / hl.sum(mt.AD),
            hl.null(hl.tfloat)
        )
    )
    # Define low-level heteroplasmy window used for NUMT warning logic.

    het_expr  = (
        (mt.allele_fraction >= 0.01) &
        (mt.allele_fraction <= 0.05)
    )
    # common_low_heteroplasmy is used as the NUMT proxy at entry level.
    numt_expr = het_expr & mt.common_low_heteroplasmy

    # Keep per-sample counters and risk tier in column annotations for metadata export.
    mt = mt.annotate_cols(
        num_hets=hl.agg.count_where(het_expr & (~mt.common_low_heteroplasmy)),
        num_numt_fp=hl.agg.count_where(numt_expr)
    )
    mt = mt.annotate_cols(
        numt_fp_risk_tier=(
            hl.case()
            .when(mt.mito_cn < 50,  "high")
            .when(mt.mito_cn < 100, "elevated")
            .default("standard")
        ),
    )

    mt = mt.annotate_cols(
        low_mtcn_numt_risk = mt.mito_cn < 100
    )

    het_expr = (
        (mt.allele_fraction >= 0.01) &
        (mt.allele_fraction <= 0.05)
    )


    # Entry-level warning: low-level heteroplasmy + common_low_heteroplasmy + low mtCN sample.
    mt = mt.annotate_entries(
        numt_fp_warning = (
            het_expr &
            mt.common_low_heteroplasmy &
            mt.low_mtcn_numt_risk
        )
    )

    # Convert row rsid set into a semicolon-delimited string for VCF export.
    mt_vcf = mt.annotate_rows(
        rsid=hl.if_else(
            hl.is_defined(mt.rsid) & (hl.len(mt.rsid) > 0),
            hl.delimit(hl.array(mt.rsid), ";"),
            hl.missing(hl.tstr)
        )
    )

    # Drop rows failing the "npg" site-level filter prior to export.
    mt_vcf = mt_vcf.filter_rows(~mt_vcf.filters.contains("npg"))
    
    # Restrict VCF row INFO payload to the fields required by downstream consumers.
    mt_vcf = mt_vcf.select_rows(
        rsid=mt_vcf.rsid,
        filters=mt_vcf.filters,
        info=hl.struct(
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
    )

    # Restrict FORMAT fields and expose NUMT warning as int32 for VCF compatibility.
    mt_vcf = mt_vcf.select_entries(
        GT=mt_vcf.GT,
        DP=mt_vcf.DP,
        AD=mt_vcf.AD,
        HL=mt_vcf.HL,
        MQ=mt_vcf.MQ,
        TLOD=mt_vcf.TLOD,
        FT=mt_vcf.FT,
        numt_fp_warning=hl.int32(mt_vcf.numt_fp_warning)
    )
    # Now, create new metadata lines for the filters that were missing from the vcf
    metadata = {
        "filter": {
            "indel_stack": {
                "Description": "Variant overlaps or is near an indel stack region."
            },
            "artifact_prone_site": {
                "Description": "Variant occurs at a site annotated as artifact-prone."
            }
        }
    }

    # export the vcf with the new additional header lines
    hl.export_vcf(
        mt_vcf,
        vcf_output_path,
        metadata=metadata,
        tabix=True
    )

    # -------------------------------------------------------------------------
    # Export sample metadata TSV
    # -------------------------------------------------------------------------
    cols_ht = mt.cols()
    # Export sample-level QC and NUMT risk fields used for reporting/traceability.
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
    sample_ht.export(sample_metadata_output_path)


if __name__ == "__main__":
    main()
