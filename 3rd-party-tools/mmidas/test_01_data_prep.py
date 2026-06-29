#!/usr/bin/env python3
"""
test_01_data_prep.py — Smoke test for 01_data_prep.py

Generates minimal synthetic Allen Brain Atlas-format CSV files,
runs the data prep script end-to-end, and asserts that the output
.h5ad has the expected structure.

Run inside the Docker image:
    python3 test_01_data_prep.py

Or locally (requires anndata, numpy, pandas, scipy, scikit-learn):
    python3 test_01_data_prep.py
"""

import os
import sys
import tempfile
import textwrap
import numpy as np
import pandas as pd
import anndata as ad

# Allow running from any directory — import the script under test directly
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)
import importlib.util

spec = importlib.util.spec_from_file_location(
    "data_prep", os.path.join(SCRIPT_DIR, "01_data_prep.py")
)
data_prep = importlib.util.module_from_spec(spec)
spec.loader.exec_module(data_prep)


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------

# A small set of gene symbols that overlap with the real selected-gene list
# (enough to test gene filtering without needing the real data files).
GENE_SYMBOLS = [
    "Malat1", "Vip", "Npy", "Ptgds", "Sst",
    "Cst3", "Apoe", "Plp1", "Bsg", "Calm2",
    "Snap25", "Aldoc", "Tmsb4x", "Penk", "Meg3",
    "Apod", "Fth1", "Tac2", "Calm1", "Crh",
    # Extra genes NOT in selected list — should be filtered out
    "FakeGene1", "FakeGene2",
]
N_GENES = len(GENE_SYMBOLS)
N_VISP = 40   # cells in VISp
N_ALM = 30    # cells in ALM
SELECTED_GENES = GENE_SYMBOLS[:20]   # first 20 are "selected"
REMOVE_CLUSTERS = ["Low Quality", "CR Lhx5", "Meis2 Adamts19"]


def _make_exon_matrix(n_cells: int, rng: np.random.Generator, sample_ids: list) -> pd.DataFrame:
    """Synthetic exon matrix: genes as rows, cells as columns.
    First column is 'gene' (gene name), rest are cell columns named by sample_id
    (must match the sample_id values in the annotation CSV).
    """
    counts = rng.integers(0, 200, size=(N_GENES, n_cells))
    df = pd.DataFrame(counts, columns=sample_ids)
    df.insert(0, "gene", GENE_SYMBOLS)
    return df


def _make_samples(n_cells: int, region: str, rng: np.random.Generator) -> pd.DataFrame:
    """Synthetic samples annotation with realistic column names."""
    classes = rng.choice(["GABAergic", "Glutamatergic", "Non-neuronal"], size=n_cells, p=[0.4, 0.4, 0.2])
    subclasses = [f"Subclass_{i % 5}" for i in range(n_cells)]

    # Assign cluster labels; sprinkle in a few that should be removed
    clusters = []
    for i in range(n_cells):
        if i % 15 == 0:
            clusters.append("Low Quality")
        elif i % 20 == 0:
            clusters.append("CR Lhx5")
        elif i % 25 == 0:
            # One of the rename cases
            clusters.append("L6b VISp Col8a1 Rprm")
        else:
            clusters.append(f"Cluster_{i % 8}")

    return pd.DataFrame({
        "sample_name": [f"{region}_sample_{i}" for i in range(n_cells)],
        "sample_id": [f"{region}_id_{i}" for i in range(n_cells)],
        "seq_batch": [f"batch_{i % 3}" for i in range(n_cells)],
        "sex": rng.choice(["M", "F"], size=n_cells),
        "brain_hemisphere": ["left"] * n_cells,
        "brain_region": [region] * n_cells,
        "brain_subregion": [f"{region}_sub"] * n_cells,
        "class": classes,
        "subclass": subclasses,
        "cluster": clusters,
        "confusion_score": rng.uniform(0, 1, size=n_cells),
    })


def _make_genes_rows() -> pd.DataFrame:
    """Full reference gene list (genes-rows.csv format)."""
    return pd.DataFrame({"gene_symbol": GENE_SYMBOLS})


def _make_selected_genes() -> pd.DataFrame:
    """Selected gene subset (genes_SS_ALM-VISp.csv format)."""
    return pd.DataFrame({"genes": SELECTED_GENES})


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

def run_smoke_test():
    rng = np.random.default_rng(seed=0)
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        # Write synthetic input CSVs
        visp_exon = os.path.join(tmpdir, "visp_exon.csv")
        visp_samples = os.path.join(tmpdir, "visp_samples.csv")
        alm_exon = os.path.join(tmpdir, "alm_exon.csv")
        alm_samples = os.path.join(tmpdir, "alm_samples.csv")
        genes_rows = os.path.join(tmpdir, "genes_rows.csv")
        selected_genes = os.path.join(tmpdir, "selected_genes.csv")
        output_h5ad = os.path.join(tmpdir, "output.h5ad")

        visp_samples_df = _make_samples(N_VISP, "VISp", rng)
        alm_samples_df = _make_samples(N_ALM, "ALM", rng)
        _make_exon_matrix(N_VISP, rng, visp_samples_df["sample_id"].tolist()).to_csv(visp_exon, index=False)
        visp_samples_df.to_csv(visp_samples, index=False)
        _make_exon_matrix(N_ALM, rng, alm_samples_df["sample_id"].tolist()).to_csv(alm_exon, index=False)
        alm_samples_df.to_csv(alm_samples, index=False)
        _make_genes_rows().to_csv(genes_rows, index=False)
        _make_selected_genes().to_csv(selected_genes, index=False)

        # Run the script via its main() function
        argv = [
            "--visp_exon_matrix", visp_exon,
            "--visp_samples",     visp_samples,
            "--alm_exon_matrix",  alm_exon,
            "--alm_samples",      alm_samples,
            "--genes_rows",       genes_rows,
            "--selected_genes",   selected_genes,
            "--output_h5ad",      output_h5ad,
            "--subsample",        "15",
        ]
        print("Running 01_data_prep.main() ...")
        try:
            data_prep.main(argv)
        except SystemExit as e:
            if e.code != 0:
                failures.append(f"main() exited with code {e.code}")

        # -----------------------------------------------------------------
        # Assertions
        # -----------------------------------------------------------------
        print("\nValidating output .h5ad ...")

        # 1. Output file exists
        if not os.path.exists(output_h5ad):
            failures.append("Output .h5ad file was not created.")
        else:
            adata = ad.read_h5ad(output_h5ad)

            # 2. Shape — n_vars must equal number of selected genes (20)
            if adata.n_vars != len(SELECTED_GENES):
                failures.append(
                    f"Expected {len(SELECTED_GENES)} vars, got {adata.n_vars}"
                )

            # 3. n_obs > 0
            if adata.n_obs == 0:
                failures.append("Output has 0 cells.")

            # 4. Subsampling respected — at most 15 cells per class × 2 classes = 30
            if adata.n_obs > 30:
                failures.append(
                    f"Subsampling to 15/class should give ≤30 cells, got {adata.n_obs}"
                )

            # 5. Required obs columns present
            required_obs = {"class", "subclass", "cluster"}
            missing_obs = required_obs - set(adata.obs.columns)
            if missing_obs:
                failures.append(f"Missing obs columns: {missing_obs}")

            # 6. No excluded clusters remain
            bad = set(REMOVE_CLUSTERS) & set(adata.obs["cluster"].unique())
            if bad:
                failures.append(f"Excluded clusters still present: {bad}")

            # 7. Only neuronal cells (no "Non-neuronal")
            non_neuronal = (adata.obs["class"] == "Non-neuronal").sum()
            if non_neuronal > 0:
                failures.append(
                    f"{non_neuronal} Non-neuronal cells remain after filtering."
                )

            # 8. Cluster rename applied
            bad_names = {"L6b VISp Col8a1 Rprm", "L6 CT ALM Nxph2 Sla"}
            found_bad = bad_names & set(adata.obs["cluster"].unique())
            if found_bad:
                failures.append(
                    f"Unrenamed cluster labels still present: {found_bad}"
                )

            # 9. var_names match selected gene list
            if list(adata.var_names) != list(SELECTED_GENES):
                failures.append(
                    f"var_names mismatch.\n"
                    f"  expected: {SELECTED_GENES[:5]}...\n"
                    f"  got:      {list(adata.var_names)[:5]}..."
                )

            # 10. X is sparse
            from scipy.sparse import issparse
            if not issparse(adata.X):
                failures.append("adata.X should be a sparse matrix.")

            print(f"  n_obs={adata.n_obs}, n_vars={adata.n_vars}")
            print(f"  obs columns: {list(adata.obs.columns)}")
            print(f"  unique clusters: {sorted(adata.obs['cluster'].unique())}")

    # -----------------------------------------------------------------
    # Report
    # -----------------------------------------------------------------
    print()
    if failures:
        print("SMOKE TEST FAILED:")
        for f in failures:
            print(f"  ✗ {f}")
        sys.exit(1)
    else:
        print("SMOKE TEST PASSED ✓")


if __name__ == "__main__":
    run_smoke_test()
