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
import os
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


def _load_exon_matrix_filtered(
    exon_matrix_path: str,
    samples_path: str,
    neuronal_classes: list,
    encoding: str = "unicode_escape",
):
    """
    Load a genes-×-cells exon-count CSV and return only the neuronal-cell
    columns, avoiding materializing the full dense matrix in memory.

    Positional alignment is used: the annotation CSV must have one row per
    cell in the same column order as the exon matrix — the standard AIBS
    format.  Column headers in the exon matrix are NOT required to match
    sample_id values (they are typically sample_name strings in AIBS data).

    The CSV layout is:
      col 0           : gene identifier (string)
      cols 1 … N      : one column per cell (in the same order as annotation rows)

    Returns
    -------
    counts : np.ndarray  shape (n_neuronal_cells, n_genes), dtype float32
    anno   : pd.DataFrame  filtered annotation rows (neuronal cells only),
             reset index.
    """
    anno = pd.read_csv(samples_path, encoding=encoding)
    neuron_mask = anno["class"].isin(neuronal_classes)
    neuronal_anno = anno[neuron_mask].reset_index(drop=True)

    # Validate positional contract: number of cell columns in the exon matrix
    # must equal number of rows in the annotation.
    with open(exon_matrix_path, "r") as fh:
        n_matrix_cols = len(fh.readline().rstrip("\n").split(",")) - 1  # subtract gene col
    if n_matrix_cols != len(anno):
        raise ValueError(
            f"Exon matrix '{os.path.basename(exon_matrix_path)}' has "
            f"{n_matrix_cols} cell columns but annotation has {len(anno)} rows. "
            f"These must be matched-pair files in identical cell order."
        )

    # Map annotation row positions → exon matrix column indices.
    # Column 0 is the gene identifier; cell at annotation row p is column p+1.
    neuronal_positions = np.where(neuron_mask.values)[0]
    keep_col_indices = [0] + [int(p) + 1 for p in neuronal_positions]

    print(
        f"  Reading {len(neuronal_positions)} / {len(anno)} "
        f"neuronal-cell columns from {os.path.basename(exon_matrix_path)} ..."
    )
    df = pd.read_csv(
        exon_matrix_path,
        usecols=keep_col_indices,
        dtype={0: str},   # gene-id column stays as string
    )

    # First column is gene identifiers; all remaining are neuronal cells.
    # shape: (n_genes, n_neuronal_cells) → transpose to (n_cells, n_genes)
    counts = df.iloc[:, 1:].values.astype(np.float32).T

    return counts, neuronal_anno


# ---------------------------------------------------------------------------
# Main logic
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    remove_clusters = [c.strip() for c in args.remove_clusters.split(",")]
    neuronal_classes = [c.strip() for c in args.neuronal_classes.split(",")]

    # ------------------------------------------------------------------
    # 1. Load raw count matrices and cell metadata
    #    Uses column-selective loading (_load_exon_matrix_filtered) to
    #    avoid materializing the full ~45K×15K dense float64 matrix
    #    (~5 GB per region). Only neuronal-cell columns are read.
    # ------------------------------------------------------------------
    print("Loading VISp data...")
    visp_counts, visp_anno = _load_exon_matrix_filtered(
        args.visp_exon_matrix, args.visp_samples, neuronal_classes
    )

    print("Loading ALM data...")
    alm_counts, alm_anno = _load_exon_matrix_filtered(
        args.alm_exon_matrix, args.alm_samples, neuronal_classes
    )

    print(
        f"Raw neuronal cell counts — VISp: {len(visp_anno)}, "
        f"ALM: {len(alm_anno)}, "
        f"total: {len(visp_anno) + len(alm_anno)}"
    )

    df_anno = pd.concat([visp_anno, alm_anno], ignore_index=True)
    total_counts = np.concatenate([visp_counts, alm_counts], axis=0)

    print(f"Combined matrix shape: {total_counts.shape}")

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
