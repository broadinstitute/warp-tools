#!/usr/bin/env python3
"""
test_03c_traversal_prep.py — smoke test for 03c_traversal_prep.py

Creates synthetic data, runs the full prerequisite chain (02b → 03a → 03c),
and asserts that all expected output files are produced with the correct keys.
"""

import json
import os
import pickle
import subprocess
import sys
import tempfile

import anndata as ad
import numpy as np
import scipy.sparse as sp


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------
N_CELLS = 60
N_GENES = 20
N_CLUSTERS = 5
SEED = 42


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
    # 03c_traversal_prep uses load_data which expects a "log1p" key
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
        h5ad_path  = os.path.join(tmpdir, "test_data.h5ad")
        out_dir    = os.path.join(tmpdir, "results")
        os.makedirs(out_dir, exist_ok=True)

        print(f"Work dir: {tmpdir}")

        # Step 1 — write synthetic .h5ad
        print("\n[1/4] Writing synthetic h5ad ...")
        _make_h5ad(h5ad_path)

        # Step 2 — train a tiny mixVAE model
        print("[2/4] Running 02b_train_mixvae.py ...")
        _run(
            [
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
            ],
            "02b_train_mixvae",
        )
        ckpt_manifest_path = os.path.join(out_dir, "checkpoints_manifest.json")
        assert os.path.exists(ckpt_manifest_path), "checkpoints_manifest.json not created"

        # Step 3 — evaluate model
        print("[3/4] Running 03a_evaluate.py ...")
        _run(
            [
                "python3", "/usr/local/03a_evaluate.py",
                "--anndata_h5ad",       h5ad_path,
                "--checkpoints_manifest", ckpt_manifest_path,
                "--output_dir",         out_dir,
                "--k_select_thr",       "0.0",
                "--batch_size",         "30",
                "--seed",               str(SEED),
            ],
            "03a_evaluate",
        )
        eval_json_path = os.path.join(out_dir, "evaluation_results.json")
        assert os.path.exists(eval_json_path), "evaluation_results.json not created"

        with open(eval_json_path) as fh:
            eval_results = json.load(fh)
        model_order = eval_results["model_order"]
        print(f"  model_order = {model_order}")

        # Step 4 — traversal prep
        print("[4/4] Running 03c_traversal_prep.py ...")
        _run(
            [
                "python3", "/usr/local/03c_traversal_prep.py",
                "--anndata_h5ad",          h5ad_path,
                "--checkpoints_manifest",  ckpt_manifest_path,
                "--evaluation_results_json", eval_json_path,
                "--output_dir",            out_dir,
                "--n_traversal_steps",     "3",
                "--seed",                  str(SEED),
            ],
            "03c_traversal_prep",
        )

        # -----------------------------------------------------------------
        # Assertions
        # -----------------------------------------------------------------
        print("\nValidating outputs ...")

        state_dir = os.path.join(out_dir, "state")
        trav_path  = os.path.join(state_dir, f"traversal_pca_K_{model_order}.pickle")
        smu_path   = os.path.join(state_dir, f"state_mu_pca_K_{model_order}.pickle")
        tax_path   = os.path.join(out_dir,   f"taxonomy_order_K_{model_order}.npy")
        mani_path  = os.path.join(out_dir,   "traversal_manifest.json")

        for p in (trav_path, smu_path, tax_path, mani_path):
            assert os.path.exists(p), f"Missing output: {p}"

        # traversal_pca pickle keys and shapes
        with open(trav_path, "rb") as fh:
            trav = pickle.load(fh)
        required_keys = {"V_g_mean", "V_g_std", "g_subset", "c_cat", "pathways"}
        assert required_keys == set(trav.keys()), (
            f"traversal pickle keys mismatch: {set(trav.keys())}"
        )
        n_cats = len(trav["c_cat"])
        assert n_cats > 0, "c_cat is empty"
        assert trav["V_g_mean"].ndim == 4, "V_g_mean should be 4-D"
        assert trav["V_g_mean"].shape[0] == n_cats, "V_g_mean axis-0 ≠ n_cats"
        assert trav["V_g_mean"].shape[-1] == N_GENES, "V_g_mean last dim ≠ n_genes"
        assert isinstance(trav["g_subset"], list), "g_subset should be a list"
        assert isinstance(trav["pathways"], list), "pathways should be a list"

        # state_mu_pca pickle keys and shapes
        with open(smu_path, "rb") as fh:
            smu = pickle.load(fh)
        required_smu_keys = {"s_mu", "s_travers", "c_cat"}
        assert required_smu_keys == set(smu.keys()), (
            f"state_mu pickle keys mismatch: {set(smu.keys())}"
        )
        assert len(smu["s_travers"]) == n_cats, "s_travers length ≠ n_cats"
        assert np.array_equal(smu["c_cat"], trav["c_cat"]), "c_cat mismatch between pickles"

        n_arm = trav["V_g_mean"].shape[1]
        n_steps = 3
        for ic in range(n_cats):
            s = smu["s_travers"][ic]
            assert s.shape == (n_arm, n_steps, 2), (
                f"s_travers[{ic}] shape {s.shape} ≠ ({n_arm}, {n_steps}, 2)"
            )

        # taxonomy order
        tax_order = np.load(tax_path)
        assert tax_order.ndim == 1, "taxonomy_order should be 1-D"
        assert len(tax_order) == n_cats, (
            f"taxonomy_order length {len(tax_order)} ≠ n_cats {n_cats}"
        )

        # manifest
        with open(mani_path) as fh:
            mani = json.load(fh)
        for key in ("model_order", "n_arm", "state_dim", "n_gene", "n_cats",
                    "traversal_pickle", "state_mu_pickle", "taxonomy_order"):
            assert key in mani, f"Missing key '{key}' in traversal_manifest.json"

        print("\nAll assertions passed. 03c_traversal_prep smoke test PASSED.")


if __name__ == "__main__":
    main()
