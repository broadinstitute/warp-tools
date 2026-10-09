"""Reproducibility check for label_transfer_from_preprocessed.py --seed.

Trains SCVI/SCANVI twice on a small synthetic query + reference with --seed 42 and asserts the
predicted labels and SCANVI embedding are identical, then runs once unseeded to confirm the default
path still works. CPU-only; CI runs it inside the built image before pushing.

    python3 test_seed_reproducibility.py [path/to/label_transfer_from_preprocessed.py]
"""
import os
import subprocess
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

script = sys.argv[1] if len(sys.argv) > 1 else "/usr/local/label_transfer_from_preprocessed.py"
rng = np.random.default_rng(0)
n_genes, n_types = 1000, 3
profiles = rng.gamma(1.0, 1.0, (n_types, n_genes))


def make(n, prefix, batch, labeled):
    types = rng.integers(n_types, size=n)
    obs = pd.DataFrame({"batch": batch}, index=[f"{prefix}{i}" for i in range(n)])
    obs["final_annotation"] = [f"type{t}" for t in types] if labeled else "Unknown"
    X = sp.csr_matrix(rng.poisson(profiles[types] * 2).astype(np.float32))
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=[f"g{j}" for j in range(n_genes)]))


def run(tmp, name, *extra):
    workdir = os.path.join(tmp, name)
    os.mkdir(workdir)
    subprocess.run([sys.executable, script, "--gex", "../gex.h5ad", "--ref", "../ref.h5ad",
                    "--input-id", "t", "--max-epochs", "2", *extra], cwd=workdir, check=True)
    return ad.read_h5ad(os.path.join(workdir, "t_SCANVI_predictions.h5ad"))


with tempfile.TemporaryDirectory() as tmp:
    make(300, "q", "query", False).write_h5ad(os.path.join(tmp, "gex.h5ad"))
    make(300, "r", "ref", True).write_h5ad(os.path.join(tmp, "ref.h5ad"))
    first = run(tmp, "seeded1", "--seed", "42")
    second = run(tmp, "seeded2", "--seed", "42")
    run(tmp, "unseeded")

same = (first.obs["celltype"].astype(str).equals(second.obs["celltype"].astype(str))
        and np.array_equal(first.obsm["X_scANVI"], second.obsm["X_scANVI"]))
print(("PASSED" if same else "FAILED") + ": --seed 42 labels and embedding "
      + ("identical" if same else "differ") + " across runs")
sys.exit(0 if same else 1)
