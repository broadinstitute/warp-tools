#!/usr/bin/env python3
"""
test_05_state_traversal.py — smoke test for 05_state_traversal.py

Runs the full upstream chain (02b → 03a → 03c → 05) with synthetic data
and asserts that all expected state-scatter figures are produced.
KEGG is omitted so only Plot A (state-space scatter) is exercised.
"""

import json
import os
import subprocess
import sys
import tempfile

import anndata as ad
import numpy as np
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------
N_CELLS = 60
N_GENES = 20
N_CLUSTERS = 4
SEED = 42
N_SEL_CATS = 3   # request only 3 categories to keep test fast


def _make_h5ad(path: str) -> None:
    rng = np.random.default_rng(SEED)
    X = rng.poisson(lam=5.0, size=(N_CELLS, N_GENES)).astype(np.float32)
    cluster_labels = rng.choice(
        [f"cluster_{i}" for i in range(N_CLUSTERS)], size=N_CELLS
    )
    adata = ad.AnnData(X=sp.csr_matrix(X))
    adata.var_names = [f"gene_{i:03d}" for i in range(N_GENES)]
    adata.obs["cluster"]  = cluster_labels
    adata.obs["subclass"] = cluster_labels
    adata.obs["class_"]   = cluster_labels
    adata.layers["log1p"] = np.log1p(X)
    adata.write_h5ad(path)


def _run(cmd: list, label: str) -> None:
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"\n=== {label} FAILED ===")
        print("STDOUT:", result.stdout[-3000:])
        print("STDERR:", result.stderr[-3000:])
        sys.exit(1)
    print(f"  {label}: OK")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path = os.path.join(tmpdir, "test_data.h5ad")
        out_dir   = os.path.join(tmpdir, "results")
        os.makedirs(out_dir, exist_ok=True)

        print(f"Work dir: {tmpdir}")
        print("\n[1/5] Writing synthetic h5ad ...")
        _make_h5ad(h5ad_path)

        print("[2/5] Running 02b_train_mixvae.py ...")
        _run([
            "python3", "/usr/local/02b_train_mixvae.py",
            "--anndata_h5ad", h5ad_path,
            "--output_dir",   out_dir,
            "--n_categories", "8",
            "--state_dim",    "2",
            "--n_arm",        "2",
            "--latent_dim",   "8",
            "--fc_dim",       "16",
            "--n_epoch",      "2",
            "--n_epoch_p",    "1",
            "--min_con",      "0.0",
            "--max_prun_it",  "1",
            "--batch_size",   "30",
            "--seed",         str(SEED),
        ], "02b_train_mixvae")
        ckpt_mf   = os.path.join(out_dir, "checkpoints_manifest.json")

        print("[3/5] Running 03a_evaluate.py ...")
        _run([
            "python3", "/usr/local/03a_evaluate.py",
            "--anndata_h5ad",         h5ad_path,
            "--checkpoints_manifest", ckpt_mf,
            "--output_dir",           out_dir,
            "--k_select_thr",         "0.0",
            "--batch_size",           "30",
            "--seed",                 str(SEED),
        ], "03a_evaluate")
        eval_json = os.path.join(out_dir, "evaluation_results.json")
        with open(eval_json) as fh:
            model_order = json.load(fh)["model_order"]
        print(f"  model_order = {model_order}")

        print("[4/5] Running 03c_traversal_prep.py (no KEGG) ...")
        _run([
            "python3", "/usr/local/03c_traversal_prep.py",
            "--anndata_h5ad",           h5ad_path,
            "--checkpoints_manifest",   ckpt_mf,
            "--evaluation_results_json", eval_json,
            "--output_dir",             out_dir,
            "--n_traversal_steps",      "3",
            "--seed",                   str(SEED),
        ], "03c_traversal_prep")
        trav_mf = os.path.join(out_dir, "traversal_manifest.json")
        with open(trav_mf) as fh:
            n_cats = json.load(fh)["n_cats"]
        n_sel = min(N_SEL_CATS, n_cats)
        print(f"  n_cats = {n_cats}, using n_selected_cats = {n_sel}")

        print("[5/5] Running 05_state_traversal.py ...")
        _run([
            "python3", "/usr/local/05_state_traversal.py",
            "--anndata_h5ad",           h5ad_path,
            "--checkpoints_manifest",   ckpt_mf,
            "--evaluation_results_json", eval_json,
            "--traversal_manifest",     trav_mf,
            "--arm",                    "0",
            "--n_selected_cats",        str(n_sel),
            "--batch_size",             "30",
            "--seed",                   str(SEED),
        ], "05_state_traversal")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        mani_path = os.path.join(out_dir, "state_traversal_manifest.json")
        assert os.path.exists(mani_path), "state_traversal_manifest.json not found"

        with open(mani_path) as fh:
            mani = json.load(fh)

        # Scatter figures.
        #
        # n_selected_cats is capped at the number of *populated* categories, not
        # just at --n_selected_cats: selecting categories with no cells assigned
        # produces figures that are all the same plot under different titles.
        # On this synthetic 2-epoch model most categories are empty, so
        # n_selected_cats is legitimately below the requested n_sel.
        assert mani["n_selected_cats"] <= n_sel, (
            f"n_selected_cats={mani['n_selected_cats']} exceeds requested {n_sel}"
        )
        assert mani["n_selected_cats"] == len(mani["selected_c"]), \
            "n_selected_cats does not match len(selected_c)"
        assert mani["n_selected_cats"] > 0, "no categories selected"

        # Every selected category must have cells, listed largest-first.
        counts = mani["selected_c_n_cells"]
        assert len(counts) == len(mani["selected_c"]), \
            "selected_c_n_cells length mismatch"
        assert all(n > 0 for n in counts), \
            f"selected categories with no cells assigned: {counts}"
        assert counts == sorted(counts, reverse=True), \
            f"selected categories not ranked largest-first: {counts}"
        assert 0 < mani["n_populated_categories"] <= mani["n_active_categories"], (
            f"n_populated_categories={mani['n_populated_categories']} outside "
            f"(0, n_active_categories={mani['n_active_categories']}]"
        )

        s_mu_dir = os.path.join(out_dir, "state", "s_mu")
        for cc in mani["selected_c"]:
            p = os.path.join(s_mu_dir, f"state_mu_arm_0_c_{cc}.png")
            assert os.path.exists(p) and os.path.getsize(p) > 0, \
                f"Missing/empty scatter: {p}"

        # No KEGG → no pathway/summary figs
        assert mani["n_pathways"] == 0, "Expected 0 pathways (no KEGG supplied)"
        assert mani["pathway_figs"] == [], "Expected empty pathway_figs"
        assert mani["summary_figs"] == [], "Expected empty summary_figs"

        print(f"\nAll {mani['n_selected_cats']} scatter figures verified "
              f"({mani['n_populated_categories']} of "
              f"{mani['n_active_categories']} categories populated).")
        print("05_state_traversal smoke test PASSED.")


if __name__ == "__main__":
    main()
