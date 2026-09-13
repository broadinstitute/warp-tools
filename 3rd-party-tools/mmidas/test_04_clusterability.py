#!/usr/bin/env python3
"""
test_04_clusterability.py — smoke test for 04_clusterability.py

Runs the full upstream chain (02b → 03a → 03b → 04) with synthetic data
and asserts that all expected figure files are produced.
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
            "--n_categories", "6",
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
        ckpt_mf = os.path.join(out_dir, "checkpoints_manifest.json")

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

        print("[4/5] Running 03b_classify.py ...")
        _run([
            "python3", "/usr/local/03b_classify.py",
            "--anndata_h5ad",           h5ad_path,
            "--checkpoints_manifest",   ckpt_mf,
            "--evaluation_results_json", eval_json,
            "--output_dir",             out_dir,
            "--n_pca",                  "8",
            "--k_fold",                 "3",
            "--batch_size",             "30",
            "--seed",                   str(SEED),
        ], "03b_classify")
        classify_mf = os.path.join(out_dir, "classify_manifest.json")
        with open(classify_mf) as fh:
            clf = json.load(fh)
        n_arm    = clf["n_arm"]
        date_tag = clf["date_tag"]

        print("[5/5] Running 04_clusterability.py ...")
        _run([
            "python3", "/usr/local/04_clusterability.py",
            "--classify_manifest", classify_mf,
            "--output_dir",        out_dir,
        ], "04_clusterability")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        expected = [
            os.path.join(out_dir, f"classAcc_RF_K_{model_order}.png"),
            os.path.join(out_dir, f"SC_K_{model_order}_{date_tag}.png"),
            os.path.join(out_dir, "conf_Ttype_pc.png"),
        ]
        for arm in range(n_arm):
            expected += [
                os.path.join(out_dir, f"conf_Ttype_lowD_arm_{arm}.png"),
                os.path.join(out_dir, f"conf_ConsType_pc_arm_{arm}.png"),
                os.path.join(out_dir, f"conf_ConsType_lowD_arm_{arm}.png"),
            ]

        for p in expected:
            assert os.path.exists(p) and os.path.getsize(p) > 0, \
                f"Missing or empty figure: {p}"

        # Check manifest
        clust_mf = os.path.join(out_dir, "clusterability_manifest.json")
        assert os.path.exists(clust_mf), "clusterability_manifest.json not found"
        with open(clust_mf) as fh:
            cm = json.load(fh)
        assert cm["model_order"] == model_order
        assert cm["n_arm"] == n_arm
        assert len(cm["figures"]) == len(expected)

        print(f"\nAll {len(expected)} figures verified.")
        print("04_clusterability smoke test PASSED.")


if __name__ == "__main__":
    main()
