"""CI check: SCP-ready outputs, run provenance, query mapping, and the input check.

Builds a small synthetic reference and an Optimus-like sample on CPU, trains a reference-only model
the way warp-scanvi-training does (make_reference_only_model.py), labels the sample with
scanvi_sample.py, and checks the outputs. Then checks that validate_inputs.py rejects a missing
file and a duplicate id. CI runs it inside the built image before pushing.

    python3 test_scp_and_query_mapping.py
"""
import os
import subprocess
import sys
import tempfile

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from make_reference_only_model import train_reference_only_model

HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(0)
genes = [f"g{j}" for j in range(500)]
profiles = rng.gamma(1.0, 1.0, (3, len(genes)))
failures = []


def synthetic_counts(n):
    types = rng.integers(3, size=n)
    return types, sp.csr_matrix(rng.poisson(profiles[types] * 3).astype(np.float32))


def check(ok, what):
    print(("ok    " if ok else "FAIL  ") + what, flush=True)
    if not ok:
        failures.append(what)


with tempfile.TemporaryDirectory() as tmp:
    os.chdir(tmp)
    types, X = synthetic_counts(400)
    ref = ad.AnnData(X, var=pd.DataFrame(index=genes), obs=pd.DataFrame(
        {"donor_id": rng.choice(["d1", "d2"], 400), "subclass": [f"type{t}" for t in types]},
        index=[f"r{i}" for i in range(400)]))
    train_reference_only_model(ref, "donor_id", "subclass", epochs=2, seed=0, out_tar="ref_only.tar.gz")

    _, X = synthetic_counts(300)
    ad.AnnData(X, var=pd.DataFrame(index=genes), obs=pd.DataFrame(
        {"star_IsCell": True}, index=[f"BC{i}" for i in range(300)])).write_h5ad("sample.h5ad")
    subprocess.run([sys.executable, os.path.join(HERE, "scanvi_sample.py"), "--input-id", "s1",
                    "--gex", os.path.abspath("sample.h5ad"), "--scanvi-model", "ref_only.tar.gz",
                    "--query-max-epochs", "2", "--seed", "0"], check=True)

    gex = ad.read_h5ad("s1_gex_annotated_matrix.h5ad")
    predictions = ad.read_h5ad("s1_SCANVI_predictions.h5ad")
    x = gex.X[:50].toarray()
    raw = gex.raw.X[:50].toarray() if gex.raw is not None else None
    check(raw is not None and np.allclose(raw, np.round(raw)), ".raw holds integer counts")
    check(not np.allclose(x, np.round(x)), ".X is normalized, not counts")
    check("counts" not in gex.layers, "no duplicate counts layer")
    check(gex.uns.get("log1p") == {}, "log1p marker has no null base (older anndata can read it)")
    check(gex.obsm["X_umap"].shape == (gex.n_obs, 2), "X_umap has 2 columns, one row per cell")
    check(gex.obsm["X_scANVI"].shape == (gex.n_obs, 30), "X_scANVI is the 30-dimensional latent")
    check(gex.uns.get("scanvi_input_id") == "s1", "uns records the input id")
    check(gex.uns.get("scanvi_input_file") == "sample.h5ad", "uns records the input file name")
    check(gex.uns.get("scanvi_model") == "ref_only.tar.gz", "uns records the model file name")
    check(set(gex.obs["final_annotation"].astype(str)) <= {"type0", "type1", "type2"},
          "labels come from the reference-only model")
    check(predictions.n_obs == gex.n_obs, "with a supplied model the predictions hold the sample only")
    check(predictions.uns.get("scanvi_input_id") == "s1", "predictions record the input id")

    result = subprocess.run([sys.executable, os.path.join(HERE, "validate_inputs.py"),
                             "--input-ids", "a", "a", "--gex", "sample.h5ad", "missing.h5ad"],
                            capture_output=True, text=True)
    check(result.returncode == 1 and "missing.h5ad" in result.stderr and "unique" in result.stderr,
          "validate_inputs.py reports a missing file and a duplicate id")

print("PASSED" if not failures else f"FAILED: {failures}", flush=True)
sys.exit(1 if failures else 0)
