#!/usr/bin/env python3
"""
01_data_prep.py — MMIDAS Data Preparation (Smart-seq Mouse ALM/VISp)

Converts raw Allen Brain Atlas Smart-seq exon count CSVs into a
normalized, filtered AnnData (.h5ad) file ready for MMIDAS training.

Steps:
  1. Load raw exon count matrices and cell metadata for VISp and ALM.
  2. Retain only neuronal cells (GABAergic and Glutamatergic classes).
  3. Concatenate both regions into a single matrix.
  4. Normalize using log CPM: log1p(counts / sum(counts) * 1e6).
  5. Filter to the selected gene set provided in --selected_genes_csv.
  6. Remove low-quality cells and rare subclasses.
  7. Rename two cell types to match the reference taxonomy.
  8. Write the result as an AnnData .h5ad file.

Usage:
  python3 01_data_prep.py \
    --visp_exon_matrix  mouse_VISp_2018-06-14_exon-matrix.csv \
    --visp_samples      mouse_VISp_2018-06-14_samples-columns.csv \
    --alm_exon_matrix   mouse_ALM_2018-06-14_exon-matrix.csv \
    --alm_samples       mouse_ALM_2018-06-14_samples-columns.csv \
    --genes_rows        mouse_ALM_2018-06-14_genes-rows.csv \
    --selected_genes    genes_SS_ALM-VISp.csv \
    --output_h5ad       Mouse_ALM-VISp_cpm.h5ad

Optional:
  --subsample N         Retain at most N cells per class (for CI smoke tests).
  --remove_clusters     Comma-separated cluster names to exclude
                        (default: "Low Quality,CR Lhx5,Meis2 Adamts19").
"""

import argparse
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from sklearn.preprocessing import normalize


# ---------------------------------------------------------------------------
# Normalization helper (inlined to avoid runtime dependency on mmidas package)
# ---------------------------------------------------------------------------

def logcpm(counts: np.ndarray, scaler: float = 1e6) -> np.ndarray:
    """Log CPM normalization.

    Parameters
    ----------
    counts:
        Cell × gene raw count matrix (dense, float or int).
    scaler:
        Counts-per-X scaler. Default 1e6 (CPM).

    Returns
    -------
    np.ndarray
        log1p-transformed, library-size-normalized matrix.
    """
    normed = normalize(counts.astype(float), axis=1, norm="l1")
    return np.log1p(normed * scaler)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Prepare Mouse Smart-seq ALM/VISp data for MMIDAS training."
    )

    # Required inputs
    parser.add_argument(
        "--visp_exon_matrix", required=True,
        help="Path to mouse_VISp_2018-06-14_exon-matrix.csv"
    )
    parser.add_argument(
        "--visp_samples", required=True,
        help="Path to mouse_VISp_2018-06-14_samples-columns.csv"
    )
    parser.add_argument(
        "--alm_exon_matrix", required=True,
        help="Path to mouse_ALM_2018-06-14_exon-matrix.csv"
    )
    parser.add_argument(
        "--alm_samples", required=True,
        help="Path to mouse_ALM_2018-06-14_samples-columns.csv"
    )
    parser.add_argument(
        "--genes_rows", required=True,
        help="Path to mouse_ALM_2018-06-14_genes-rows.csv (full gene list)"
    )
    parser.add_argument(
        "--selected_genes", required=True,
        help="Path to genes_SS_ALM-VISp.csv (selected gene subset)"
    )
    parser.add_argument(
        "--output_h5ad", required=True,
        help="Path for the output .h5ad file"
    )

    # Optional filters
    parser.add_argument(
        "--remove_clusters",
        default="Low Quality,CR Lhx5,Meis2 Adamts19",
        help=(
            "Comma-separated list of cluster names to exclude. "
            "Default: 'Low Quality,CR Lhx5,Meis2 Adamts19'"
        )
    )
    parser.add_argument(
        "--neuronal_classes",
        default="GABAergic,Glutamatergic",
        help=(
            "Comma-separated list of cell classes to retain. "
            "Default: 'GABAergic,Glutamatergic'"
        )
    )

    # CI / smoke-test option
    parser.add_argument(
        "--subsample", type=int, default=0,
        help=(
            "If > 0, subsample to at most this many cells per neuronal class "
            "after filtering. Useful for fast CI smoke tests."
        )
    )

    return parser.parse_args(argv)


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    remove_clusters = [c.strip() for c in args.remove_clusters.split(",")]
    neuronal_classes = [c.strip() for c in args.neuronal_classes.split(",")]

    # ------------------------------------------------------------------
    # 1. Load raw count matrices and cell metadata
    # ------------------------------------------------------------------
    print("Loading VISp data...")
    df_visp_exon = pd.read_csv(args.visp_exon_matrix)
    df_visp_anno = pd.read_csv(args.visp_samples, encoding="unicode_escape")

    print("Loading ALM data...")
    df_alm_exon = pd.read_csv(args.alm_exon_matrix)
    df_alm_anno = pd.read_csv(args.alm_samples, encoding="unicode_escape")

    print(
        f"Raw cell counts — VISp: {len(df_visp_anno)}, ALM: {len(df_alm_anno)}"
    )

    # ------------------------------------------------------------------
    # 2. Filter to neuronal cells
    #    The exon matrix CSVs have genes as rows and cells as columns;
    #    the first column is the gene identifier, so we slice [:, 1:].
    # ------------------------------------------------------------------
    visp_neuron_mask = df_visp_anno["class"].isin(neuronal_classes)
    alm_neuron_mask = df_alm_anno["class"].isin(neuronal_classes)

    # shape: (n_cells, n_genes) after transpose
    visp_counts = df_visp_exon.values[:, 1:][:, visp_neuron_mask].T
    alm_counts = df_alm_exon.values[:, 1:][:, alm_neuron_mask].T

    df_anno = pd.concat(
        [df_visp_anno[visp_neuron_mask], df_alm_anno[alm_neuron_mask]],
        ignore_index=True,
    )
    total_counts = np.concatenate([visp_counts, alm_counts], axis=0)

    print(
        f"Neuronal cells retained — VISp: {visp_neuron_mask.sum()}, "
        f"ALM: {alm_neuron_mask.sum()}, "
        f"total: {total_counts.shape[0]}"
    )

    # ------------------------------------------------------------------
    # 3. Log CPM normalization
    # ------------------------------------------------------------------
    print("Normalizing to log CPM...")
    log_cpm = logcpm(total_counts)
    print(f"  Row sums after normalization (should be ~log(1e6+1)≈13.8): "
          f"min={log_cpm.sum(axis=1).min():.2f}, "
          f"max={log_cpm.sum(axis=1).max():.2f}")

    # ------------------------------------------------------------------
    # 4. Filter to selected gene set
    # ------------------------------------------------------------------
    print("Loading gene lists...")
    ref_genes_df = pd.read_csv(args.genes_rows)
    slc_genes_df = pd.read_csv(args.selected_genes)
    selected_genes = slc_genes_df["genes"].values

    # Build index mapping: for each selected gene, find its position in the
    # full gene list (ref_genes_df.gene_symbol).
    ref_symbols = ref_genes_df["gene_symbol"].values
    gene_indices = []
    missing = []
    for gene in selected_genes:
        hits = np.where(ref_symbols == gene)[0]
        if len(hits) == 0:
            missing.append(gene)
        else:
            gene_indices.append(hits[0])

    if missing:
        print(
            f"WARNING: {len(missing)} selected genes not found in reference "
            f"gene list and will be skipped: {missing[:10]}{'...' if len(missing) > 10 else ''}"
        )

    gene_indices = np.array(gene_indices, dtype=int)
    log_cpm_filtered = log_cpm[:, gene_indices]
    retained_genes = selected_genes[[i for i, g in enumerate(selected_genes) if g not in missing]]

    print(
        f"Gene filtering: {len(ref_symbols)} total → "
        f"{len(retained_genes)} selected genes retained"
    )

    # ------------------------------------------------------------------
    # 5. Remove low-quality cells and excluded cluster types
    # ------------------------------------------------------------------
    keep_mask = ~df_anno["cluster"].isin(remove_clusters)
    df_anno_filtered = df_anno[keep_mask].reset_index(drop=True)
    log_cpm_filtered = log_cpm_filtered[keep_mask.values, :]

    print(
        f"Cluster filtering: removed {(~keep_mask).sum()} cells "
        f"({remove_clusters}). Remaining: {log_cpm_filtered.shape[0]} cells."
    )

    # ------------------------------------------------------------------
    # 6. Optional subsampling (for CI smoke tests)
    # ------------------------------------------------------------------
    if args.subsample > 0:
        print(f"Subsampling to at most {args.subsample} cells per class...")
        keep_indices = []
        for cls in neuronal_classes:
            cls_idx = np.where(df_anno_filtered["class"].values == cls)[0]
            if len(cls_idx) > args.subsample:
                rng = np.random.default_rng(seed=42)
                cls_idx = rng.choice(cls_idx, size=args.subsample, replace=False)
            keep_indices.append(cls_idx)
        keep_indices = np.sort(np.concatenate(keep_indices))
        df_anno_filtered = df_anno_filtered.iloc[keep_indices].reset_index(drop=True)
        log_cpm_filtered = log_cpm_filtered[keep_indices, :]
        print(f"  After subsampling: {log_cpm_filtered.shape[0]} cells.")

    # ------------------------------------------------------------------
    # 7. Rename two cell types to match the reference taxonomy
    # ------------------------------------------------------------------
    df_anno_filtered["cluster"] = df_anno_filtered["cluster"].replace({
        "L6b VISp Col8a1 Rprm": "L6b Col8a1 Rprm",
        "L6 CT ALM Nxph2 Sla": "L6 CT Nxph2 Sla",
    })

    # ------------------------------------------------------------------
    # 8. Build AnnData and write .h5ad
    # ------------------------------------------------------------------
    print("Building AnnData object...")
    obs_cols = [
        "sample_name", "sample_id", "seq_batch", "sex",
        "brain_hemisphere", "brain_region", "brain_subregion",
        "class", "subclass", "cluster", "confusion_score",
    ]
    # Keep only columns that actually exist in the annotation dataframe
    obs_cols = [c for c in obs_cols if c in df_anno_filtered.columns]
    obs_df = df_anno_filtered[obs_cols]

    adata = ad.AnnData(X=csr_matrix(log_cpm_filtered), obs=obs_df)
    adata.var_names = retained_genes
    if "sample_id" in obs_df.columns:
        adata.obs_names = obs_df["sample_id"].values
    else:
        adata.obs_names = [str(i) for i in range(len(obs_df))]

    print(
        f"Final AnnData shape: {adata.n_obs} cells × {adata.n_vars} genes"
    )
    print(f"  obs columns: {list(adata.obs.columns)}")
    print(f"  Unique clusters: {adata.obs['cluster'].nunique()}")

    print(f"Writing output to: {args.output_h5ad}")
    adata.write_h5ad(args.output_h5ad)
    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
