#!/usr/bin/env python3
"""
test_03a_evaluate.py — Smoke test for 03a_evaluate.py

Runs the full 02b training (2 epochs, 1 prune) on a tiny synthetic .h5ad,
then runs 03a_evaluate on the resulting checkpoints and asserts that all
expected outputs (JSON, pickle, figures) are produced.

Run inside the Docker image:
    python3 /usr/local/test_03a_evaluate.py
"""

import importlib.util
import json
import os
import pickle
import re
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

def _load_script(name, filename):
    spec = importlib.util.spec_from_file_location(
        name, os.path.join(SCRIPT_DIR, filename)
    )
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

train_mixvae = _load_script("train_mixvae", "02b_train_mixvae.py")
evaluate     = _load_script("evaluate",     "03a_evaluate.py")


# ---------------------------------------------------------------------------
# Synthetic data
# ---------------------------------------------------------------------------

N_CELLS = 80
N_GENES = 20
N_CATEGORIES = 5   # keep small so pruning is fast


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
# Regression test: model_order -> checkpoint resolution
#
# Guards two bugs that together caused evaluation_results.json to report a
# model_order taken from a checkpoint nothing had selected:
#   1. model_order == n_categories yields pruning round 0, but the round-0 model
#      is named "before_pruning" — the glob for "after_pruning_0_*" matched
#      nothing and control fell through to a fallback.
#   2. That fallback used sorted(all_after)[-1], a lexicographic sort in which
#      round 9 sorts after round 14.
# ---------------------------------------------------------------------------

CKPT_TS = "2026-01-01-00-00-00"


def _retained(path, n_categories):
    if "before_pruning" in os.path.basename(path):
        return n_categories
    return n_categories - int(
        re.search(r"after_pruning_(\d+)_", os.path.basename(path)).group(1)
    )


def run_checkpoint_resolution_test():
    print("Checkpoint resolution regression test ...")
    failures = []
    n_categories, max_prun_it = 120, 14

    with tempfile.TemporaryDirectory() as tmpdir:
        model_dir = os.path.join(tmpdir, "model")
        os.makedirs(model_dir)
        manifest_path = os.path.join(tmpdir, "checkpoints_manifest.json")
        with open(manifest_path, "w") as fh:
            json.dump({"n_categories": n_categories}, fh)

        open(os.path.join(model_dir,
             f"cpl_mixVAE_model_before_pruning_{CKPT_TS}.pth"), "w").close()
        for r in range(1, max_prun_it + 1):
            open(os.path.join(model_dir,
                 f"cpl_mixVAE_model_after_pruning_{r}_{CKPT_TS}.pth"), "w").close()

        # Every reachable model_order must resolve to a checkpoint retaining
        # exactly that many categories — including the unpruned endpoint, which
        # is the case that previously silently returned round 9.
        for order in range(n_categories - max_prun_it, n_categories + 1):
            sel = evaluate.resolve_checkpoint(manifest_path, n_categories, order)
            if sel is None:
                failures.append(f"model_order={order}: no checkpoint resolved")
                continue
            got = _retained(sel, n_categories)
            if got != order:
                failures.append(
                    f"model_order={order}: resolved checkpoint retains {got} "
                    f"categories ({os.path.basename(sel)})"
                )
        print(f"  checked model_order "
              f"{n_categories - max_prun_it}..{n_categories}")

        # An unreachable model_order must fall back to the deepest available
        # pruning round rather than an arbitrary one.
        sel = evaluate.resolve_checkpoint(manifest_path, n_categories, 50)
        if sel is None or _retained(sel, n_categories) != n_categories - max_prun_it:
            failures.append(
                f"unreachable model_order=50 resolved to "
                f"{os.path.basename(sel) if sel else None}, expected pruning "
                f"round {max_prun_it}"
            )

        # No checkpoints at all: resolve, don't crash.
        empty_dir = os.path.join(tmpdir, "empty")
        os.makedirs(os.path.join(empty_dir, "model"))
        empty_manifest = os.path.join(empty_dir, "checkpoints_manifest.json")
        with open(empty_manifest, "w") as fh:
            json.dump({"n_categories": n_categories}, fh)
        if evaluate.resolve_checkpoint(empty_manifest, n_categories, 111) is not None:
            failures.append("empty model dir should resolve to None")

    return failures


# ---------------------------------------------------------------------------
# Smoke test
# ---------------------------------------------------------------------------

def run_smoke_test():
    rng = np.random.default_rng(seed=99)
    failures = run_checkpoint_resolution_test()

    with tempfile.TemporaryDirectory() as tmpdir:
        h5ad_path  = os.path.join(tmpdir, "synthetic.h5ad")
        train_dir  = os.path.join(tmpdir, "run_train")
        eval_dir   = os.path.join(tmpdir, "run_eval")

        _make_h5ad(h5ad_path, rng)
        print(f"Synthetic .h5ad: {N_CELLS} cells × {N_GENES} genes")

        # ------------------------------------------------------------------
        # Step 1: train (reuse 02b smoke-test settings)
        # ------------------------------------------------------------------
        train_argv = [
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
        ]
        print("\nRunning 02b_train_mixvae.main() ...")
        try:
            train_mixvae.main(train_argv)
        except SystemExit as e:
            if e.code not in (None, 0):
                failures.append(f"train main() exited with code {e.code}")

        manifest_path = os.path.join(train_dir, "checkpoints_manifest.json")
        if not os.path.exists(manifest_path):
            failures.append("checkpoints_manifest.json not produced by training — cannot proceed.")
            _report(failures)

        # ------------------------------------------------------------------
        # Step 2: evaluate
        # ------------------------------------------------------------------
        eval_argv = [
            "--anndata_h5ad",         h5ad_path,
            "--checkpoints_manifest", manifest_path,
            "--output_dir",           eval_dir,
            "--k_select_thr",         "0.0",   # permissive: any consensus passes
            "--batch_size",           "16",
            "--seed",                 "0",
        ]
        print("\nRunning 03a_evaluate.main() ...")
        try:
            evaluate.main(eval_argv)
        except SystemExit as e:
            if e.code not in (None, 0):
                failures.append(f"evaluate main() exited with code {e.code}")

        # ------------------------------------------------------------------
        # Assertions
        # ------------------------------------------------------------------
        print("\nValidating outputs ...")

        results_path = os.path.join(eval_dir, "evaluation_results.json")

        # 1. evaluation_results.json exists
        if not os.path.exists(results_path):
            failures.append(f"evaluation_results.json not found: {results_path}")
        else:
            with open(results_path) as fh:
                results = json.load(fh)

            # 2. Required keys, including the review fields a human needs to
            #    tell an accepted model from a fallback (see resolve_checkpoint).
            for key in ("model_order", "selected_model", "summary_pickle",
                        "avg_consensus", "figures",
                        "k_selection_met_threshold",
                        "k_selection_suggested_model_order",
                        "n_populated_categories",
                        "n_populated_categories_per_arm",
                        "collapse_warning"):
                if key not in results:
                    failures.append(f"Results JSON missing key: '{key}'")

            # 2b. n_populated_categories must never exceed model_order, and the
            #     per-arm list must have one entry per arm.
            n_pop = results.get("n_populated_categories", -1)
            per_arm = results.get("n_populated_categories_per_arm", [])
            if not (0 < n_pop <= results.get("model_order", 0)):
                failures.append(
                    f"n_populated_categories={n_pop} outside "
                    f"(0, model_order={results.get('model_order')}]"
                )
            if len(per_arm) != 2:
                failures.append(
                    f"n_populated_categories_per_arm has {len(per_arm)} "
                    f"entries, expected 2 (n_arm)"
                )
            else:
                print(f"  n_populated_categories = {n_pop} (per arm: {per_arm})")

            # 3. model_order is a positive integer <= n_categories
            mo = results.get("model_order", 0)
            if not (0 < mo <= N_CATEGORIES):
                failures.append(
                    f"model_order={mo} outside valid range (0, {N_CATEGORIES}]."
                )
            else:
                print(f"  model_order = {mo}")

            # 4. selected_model file exists
            sel = results.get("selected_model", "")
            if not os.path.exists(sel):
                failures.append(f"selected_model path does not exist: {sel}")

            # 5. summary pickle exists and is loadable
            sp = results.get("summary_pickle", "")
            if not os.path.exists(sp):
                failures.append(f"summary_pickle not found: {sp}")
            else:
                with open(sp, "rb") as fh:
                    summary = pickle.load(fh)
                for key in ("con_mean", "num_pruned"):
                    if key not in summary:
                        failures.append(f"Summary pickle missing key: '{key}'")
                print(f"  summary_pickle keys: {list(summary.keys())}")

            # 6. All listed figure files exist and are non-empty
            for fig_path in results.get("figures", []):
                if not os.path.exists(fig_path):
                    failures.append(f"Figure file not found: {fig_path}")
                elif os.path.getsize(fig_path) == 0:
                    failures.append(f"Figure file is empty: {fig_path}")
                else:
                    print(f"  figure: {os.path.basename(fig_path)}")

            # 7. avg_consensus is a float in [0, 1]
            ac = results.get("avg_consensus", -1)
            if not (0.0 <= ac <= 1.0):
                failures.append(f"avg_consensus={ac} not in [0, 1].")
            else:
                print(f"  avg_consensus = {ac:.4f}")

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
