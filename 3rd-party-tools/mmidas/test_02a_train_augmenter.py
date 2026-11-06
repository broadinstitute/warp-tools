#!/usr/bin/env python3
"""
test_02a_train_augmenter.py — Smoke test for 02a_train_augmenter.py

Creates a minimal synthetic .h5ad file, runs the augmenter training script
for 2 epochs, and validates that the expected outputs (.pth checkpoint,
loss_curve.png, augmenter_manifest.json) were produced.

Run inside the Docker image:
    python3 /usr/local/test_02a_train_augmenter.py

Or locally (requires torch, anndata, numpy, scipy, scikit-learn, mmidas):
    python3 test_02a_train_augmenter.py
"""

import importlib.util
import json
import os
import sys
import tempfile

import numpy as np
import anndata as ad
import pandas as pd
from scipy.sparse import csr_matrix

# ---------------------------------------------------------------------------
# Import script under test
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "train_augmenter", os.path.join(SCRIPT_DIR, "02a_train_augmenter.py")
)
train_augmenter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train_augmenter)


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------

N_CELLS = 64    # must be > batch_size=16 with drop_last=True (4 full batches)
N_GENES = 20


def _make_synthetic_h5ad(path: str, rng: np.random.Generator) -> None:
    """Write a minimal .h5ad with log-CPM-like expression."""
    X = rng.uniform(0, 5, size=(N_CELLS, N_GENES)).astype(np.float32)
    obs = pd.DataFrame(
        {
            "cluster":   [f"Cluster_{i % 4}" for i in range(N_CELLS)],
            "subclass":  [f"Sub_{i % 2}"     for i in range(N_CELLS)],
            "class":     ["GABAergic" if i % 2 == 0 else "Glutamatergic"
                          for i in range(N_CELLS)],
        },
        index=[f"cell_{i}" for i in range(N_CELLS)],
    )
    var = pd.DataFrame(index=[f"gene_{i}" for i in range(N_GENES)])
    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    adata.write_h5ad(path)


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def run_smoke_test():
    rng = np.random.default_rng(seed=7)
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path   = os.path.join(tmpdir, "synthetic.h5ad")
        output_dir  = os.path.join(tmpdir, "augmenter_out")

        _make_synthetic_h5ad(h5ad_path, rng)
        print(f"Synthetic .h5ad written: {h5ad_path}  ({N_CELLS} cells × {N_GENES} genes)")

        argv = [
            "--anndata_h5ad", h5ad_path,
            "--output_dir",   output_dir,
            "--tag",          "smoketest",
            # Tiny model for speed
            "--z_dim",        "4",
            "--noise_dim",    "8",
            "--fc_dim",       "16",
            "--x_drop",       "0.0",
            # Minimal training
            "--n_epoch",      "2",
            "--batch_size",   "16",
            "--lr",           "0.001",
        ]

        print("\nRunning 02a_train_augmenter.main() ...")
        try:
            train_augmenter.main(argv)
        except SystemExit as e:
            if e.code not in (None, 0):
                failures.append(f"main() exited with non-zero code: {e.code}")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        manifest_path = os.path.join(output_dir, "augmenter_manifest.json")

        # 1. output directory created
        if not os.path.isdir(output_dir):
            failures.append(f"Output directory not created: {output_dir}")

        # 2. manifest exists and is valid JSON
        if not os.path.exists(manifest_path):
            failures.append(f"augmenter_manifest.json not found at: {manifest_path}")
        else:
            with open(manifest_path) as fh:
                manifest = json.load(fh)

            # 3. required manifest keys
            for key in ("output_dir", "checkpoint", "n_gene", "z_dim", "noise_dim"):
                if key not in manifest:
                    failures.append(f"Manifest missing key: '{key}'")

            # 4. checkpoint file exists on disk
            ckpt = manifest.get("checkpoint", "")
            if not os.path.exists(ckpt):
                failures.append(f"Checkpoint file not found: {ckpt}")
            else:
                size = os.path.getsize(ckpt)
                print(f"  Checkpoint: {os.path.basename(ckpt)}  ({size:,} bytes)")
                if size == 0:
                    failures.append("Checkpoint file is empty (0 bytes).")

            # 5. checkpoint filename contains our tag
            if "smoketest" not in os.path.basename(ckpt):
                failures.append(
                    f"Expected 'smoketest' in checkpoint filename, got: {os.path.basename(ckpt)}"
                )

            # 6. n_gene matches synthetic data
            if manifest.get("n_gene") != N_GENES:
                failures.append(
                    f"Manifest n_gene={manifest.get('n_gene')}, expected {N_GENES}."
                )

        # 7. loss_curve.png produced
        loss_curve = os.path.join(output_dir, "loss_curve.png")
        if not os.path.exists(loss_curve):
            failures.append("loss_curve.png not found in output directory.")
        else:
            print(f"  loss_curve.png: {os.path.getsize(loss_curve):,} bytes")

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
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
