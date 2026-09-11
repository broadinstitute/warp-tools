#!/usr/bin/env python3
"""
test_03b_classify.py — Smoke test for 03b_classify.py

Runs the full pipeline: 02b train → 03a evaluate → 03b classify on a tiny
synthetic dataset, and validates that all expected pickle files are produced
with the correct internal structure.

Run inside the Docker image:
    python3 /usr/local/test_03b_classify.py
"""

import importlib.util
import json
import os
import pickle
import sys
import tempfile

import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

# ---------------------------------------------------------------------------
# Import scripts under test
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

def _load(name, filename):
    spec = importlib.util.spec_from_file_location(name, os.path.join(SCRIPT_DIR, filename))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

train_mixvae = _load("train_mixvae", "02b_train_mixvae.py")
evaluate     = _load("evaluate",     "03a_evaluate.py")
classify     = _load("classify",     "03b_classify.py")


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------

N_CELLS      = 80
N_GENES      = 20
N_CATEGORIES = 5


def _make_h5ad(path, rng):
    X = rng.uniform(0, 5, size=(N_CELLS, N_GENES)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cluster":  [f"Cluster_{i % N_CATEGORIES}" for i in range(N_CELLS)],
            "subclass": [f"Sub_{i % 2}"                for i in range(N_CELLS)],
            "class":    ["GABAergic" if i % 2 == 0 else "Glutamatergic"
                         for i in range(N_CELLS)],
        },
        index=[f"cell_{i}" for i in range(N_CELLS)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(N_GENES)])
    ad.AnnData(X=csr_matrix(X), obs=obs, var=var).write_h5ad(path)


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def run_smoke_test():
    rng = np.random.default_rng(seed=77)
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path = os.path.join(tmpdir, "synthetic.h5ad")
        train_dir = os.path.join(tmpdir, "run_train")
        eval_dir  = os.path.join(tmpdir, "run_eval")
        cls_dir   = os.path.join(tmpdir, "run_eval")   # classify writes into eval dir

        _make_h5ad(h5ad_path, rng)

        # ------------------------------------------------------------------
        # Step 1: train
        # ------------------------------------------------------------------
        train_mixvae.main([
            "--anndata_h5ad", h5ad_path,
            "--output_dir",   train_dir,
            "--n_categories", str(N_CATEGORIES),
            "--state_dim",    "1",
            "--n_arm",        "2",
            "--latent_dim",   "3",
            "--fc_dim",       "16",
            "--x_drop",       "0.0",
            "--tau",          "0.2",
            "--n_epoch",      "2",
            "--n_epoch_p",    "1",
            "--max_prun_it",  "1",
            "--batch_size",   "16",
            "--seed",         "0",
        ])

        manifest_path = os.path.join(train_dir, "checkpoints_manifest.json")

        # ------------------------------------------------------------------
        # Step 2: evaluate
        # ------------------------------------------------------------------
        evaluate.main([
            "--anndata_h5ad",         h5ad_path,
            "--checkpoints_manifest", manifest_path,
            "--output_dir",           eval_dir,
            "--k_select_thr",         "0.0",
            "--batch_size",           "16",
            "--seed",                 "0",
        ])

        eval_results_path = os.path.join(eval_dir, "evaluation_results.json")

        # ------------------------------------------------------------------
        # Step 3: classify
        # ------------------------------------------------------------------
        print("\nRunning 03b_classify.main() ...")
        try:
            classify.main([
                "--anndata_h5ad",            h5ad_path,
                "--checkpoints_manifest",    manifest_path,
                "--evaluation_results_json", eval_results_path,
                "--output_dir",              cls_dir,
                "--n_pca",                   "5",
                "--k_fold",                  "2",
                "--batch_size",              "16",
                "--seed",                    "0",
                "--date_tag",                "smoketest",
            ])
        except SystemExit as e:
            if e.code not in (None, 0):
                failures.append(f"classify main() exited with code {e.code}")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        cls_manifest_path = os.path.join(cls_dir, "classify_manifest.json")
        if not os.path.exists(cls_manifest_path):
            failures.append(f"classify_manifest.json not found: {cls_manifest_path}")
            _report(failures)

        with open(cls_manifest_path) as fh:
            cls_manifest = json.load(fh)

        # 1. Manifest has expected keys
        for key in ("pickles", "model_order", "n_ttype", "n_pca", "k_fold"):
            if key not in cls_manifest:
                failures.append(f"Manifest missing key: '{key}'")

        # 2. All listed pickles exist and have correct structure
        pickles = cls_manifest.get("pickles", [])
        if not pickles:
            failures.append("No pickle files listed in manifest.")

        print(f"  {len(pickles)} pickle(s) to validate:")
        for pkl_path in pickles:
            if not os.path.exists(pkl_path):
                failures.append(f"Pickle file not found: {pkl_path}")
                continue

            with open(pkl_path, "rb") as fh:
                d = pickle.load(fh)

            name = os.path.basename(pkl_path)
            print(f"    {name}")

            # 3. Required keys
            for key in ("acc_T_adj", "sc_T", "conf_mat"):
                if key not in d:
                    failures.append(f"{name}: missing key '{key}'")

            # 4. acc_T_adj is an array of length k_fold
            acc = d.get("acc_T_adj", np.array([]))
            if len(acc) != 2:  # k_fold=2
                failures.append(f"{name}: acc_T_adj length={len(acc)}, expected 2.")

            # 5. All accuracies in [0, 1]
            if len(acc) > 0 and not np.all((acc >= 0) & (acc <= 1)):
                failures.append(f"{name}: acc_T_adj contains values outside [0,1].")

            # 6. sc_T is a non-empty list
            sc = d.get("sc_T", [])
            if not isinstance(sc, list) or len(sc) == 0:
                failures.append(f"{name}: sc_T should be a non-empty list.")

            # 7. conf_mat is a 2-D square array
            cm = d.get("conf_mat", np.array([]))
            if cm.ndim != 2 or cm.shape[0] != cm.shape[1]:
                failures.append(f"{name}: conf_mat shape {cm.shape} is not square 2-D.")

    _report(failures)


def _report(failures):
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
