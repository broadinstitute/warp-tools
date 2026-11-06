#!/usr/bin/env python3
"""
test_02b_train_mixvae.py — Smoke test for 02b_train_mixvae.py

Creates a minimal synthetic .h5ad file, runs the training script for
2 epochs + 1 pruning iteration, and validates that the expected outputs
(checkpoint .pth files and checkpoints_manifest.json) were produced.

Run inside the Docker image:
    python3 /usr/local/test_02b_train_mixvae.py

Or locally (requires torch, anndata, numpy, scipy, scikit-learn, mmidas):
    python3 test_02b_train_mixvae.py
"""

import importlib.util
import json
import os
import sys
import tempfile

import numpy as np
import anndata as ad
from scipy.sparse import csr_matrix

# ---------------------------------------------------------------------------
# Import script under test
# ---------------------------------------------------------------------------
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "train_mixvae", os.path.join(SCRIPT_DIR, "02b_train_mixvae.py")
)
train_mixvae = importlib.util.module_from_spec(spec)
spec.loader.exec_module(train_mixvae)


# ---------------------------------------------------------------------------
# Synthetic data helpers
# ---------------------------------------------------------------------------

N_CELLS = 80      # large enough for batch_size=16 with drop_last=True (4 full batches)
N_GENES = 20
N_CLUSTERS = 4    # distinct cluster labels in obs


def _make_synthetic_h5ad(path: str, rng: np.random.Generator) -> None:
    """Write a minimal .h5ad with log-CPM-like expression and required obs columns."""
    # Synthetic log-CPM expression: all positive floats
    X = rng.uniform(0, 5, size=(N_CELLS, N_GENES)).astype(np.float32)

    clusters = np.array([f"Cluster_{i % N_CLUSTERS}" for i in range(N_CELLS)])
    subclasses = np.array([f"Subclass_{i % 2}" for i in range(N_CELLS)])
    classes = np.array(["GABAergic" if i % 2 == 0 else "Glutamatergic" for i in range(N_CELLS)])

    import pandas as pd
    obs = pd.DataFrame(
        {"cluster": clusters, "subclass": subclasses, "class": classes},
        index=[f"cell_{i}" for i in range(N_CELLS)],
    )
    var = pd.DataFrame(
        index=[f"gene_{i}" for i in range(N_GENES)]
    )

    adata = ad.AnnData(X=csr_matrix(X), obs=obs, var=var)
    adata.write_h5ad(path)


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def run_smoke_test():
    rng = np.random.default_rng(seed=42)
    failures = []

    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path = os.path.join(tmpdir, "synthetic.h5ad")
        output_dir = os.path.join(tmpdir, "run_smoke")

        # Write synthetic input
        _make_synthetic_h5ad(h5ad_path, rng)
        print(f"Synthetic .h5ad written: {h5ad_path}  ({N_CELLS} cells × {N_GENES} genes)")

        # Run training with minimal settings
        argv = [
            "--anndata_h5ad", h5ad_path,
            "--output_dir",   output_dir,
            # Architecture: tiny model for speed
            "--n_categories",  "5",
            "--state_dim",     "1",
            "--n_arm",         "2",
            "--latent_dim",    "3",
            "--fc_dim",        "16",
            "--x_drop",        "0.0",   # no dropout for tiny test
            "--s_drop",        "0.0",
            "--tau",           "0.2",
            # Training: minimal epochs so the test finishes in seconds
            "--n_epoch",       "2",
            "--n_epoch_p",     "1",
            "--max_prun_it",   "1",
            "--batch_size",    "16",
            "--lr",            "0.001",
            "--seed",          "0",
        ]

        print("\nRunning 02b_train_mixvae.main() ...")
        try:
            train_mixvae.main(argv)
        except SystemExit as e:
            if e.code not in (None, 0):
                failures.append(f"main() exited with non-zero code: {e.code}")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        model_dir = os.path.join(output_dir, "model")
        manifest_path = os.path.join(output_dir, "checkpoints_manifest.json")

        # 1. model/ subdirectory was created
        if not os.path.isdir(model_dir):
            failures.append(f"model/ directory not created: {model_dir}")
        else:
            # 2. At least one .pth checkpoint exists
            pth_files = [f for f in os.listdir(model_dir) if f.endswith(".pth")]
            if not pth_files:
                failures.append("No .pth checkpoint files found in model/.")
            else:
                print(f"  Checkpoints found ({len(pth_files)}):")
                for f in sorted(pth_files):
                    fpath = os.path.join(model_dir, f)
                    print(f"    {f}  ({os.path.getsize(fpath):,} bytes)")

            # 3. before_pruning checkpoint present (n_epoch=2 > 0)
            before_pth = [f for f in pth_files if "before_pruning" in f]
            if not before_pth:
                failures.append(
                    "Expected a 'cpl_mixVAE_model_before_pruning_*.pth' file — not found."
                )

        # 4. Manifest exists and is valid JSON
        if not os.path.exists(manifest_path):
            failures.append(f"checkpoints_manifest.json not found at: {manifest_path}")
        else:
            with open(manifest_path) as fh:
                manifest = json.load(fh)

            # 5. Manifest has expected keys
            for key in ("output_dir", "checkpoints", "n_categories", "n_gene"):
                if key not in manifest:
                    failures.append(f"Manifest missing key: '{key}'")

            # 6. Manifest checkpoints list is non-empty and files exist
            if not manifest.get("checkpoints"):
                failures.append("Manifest 'checkpoints' list is empty.")
            else:
                for ckpt in manifest["checkpoints"]:
                    if not os.path.exists(ckpt):
                        failures.append(f"Manifest references missing file: {ckpt}")

            # 7. n_gene matches our synthetic data
            if manifest.get("n_gene") != N_GENES:
                failures.append(
                    f"Manifest n_gene={manifest.get('n_gene')}, expected {N_GENES}."
                )

            print(f"  Manifest: {len(manifest.get('checkpoints', []))} checkpoint(s) listed.")

        # 8. Learning curve PNG produced (cpl_mixvae saves one in model/)
        if os.path.isdir(model_dir):
            png_files = [f for f in os.listdir(model_dir) if f.endswith(".png")]
            if not png_files:
                failures.append("No learning curve .png found in model/.")
            else:
                print(f"  Learning curve figures: {png_files}")

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
