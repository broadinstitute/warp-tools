#!/usr/bin/env python3
"""
03a_evaluate.py — MMIDAS Model Evaluation (Categorical)

Loads all checkpoint .pth files produced by 02b_train_mixvae.py, runs
summarize_inference() across them, selects the optimal number of clusters
(model_order) via K_selection(), and writes figures + a JSON result.

Steps:
  1. Load the .h5ad file and build a DataLoader.
  2. Instantiate cpl_mixVAE and run summarize_inference() over all checkpoints,
     saving summary_performance_K_<K>_narm_<N>.p to the output directory.
  3. Run K_selection() to determine model_order.
  4. Re-run summarize_inference() on the selected model checkpoint alone,
     against the full dataset, to produce per-cell inference outputs.
  5. Save all figures (consensus matrix, normalised consensus, state scatter).
  6. Write evaluation_results.json containing model_order and selected_model path.

Usage:
  python3 03a_evaluate.py \\
    --anndata_h5ad       Mouse_ALM-VISp_cpm.h5ad \\
    --checkpoints_manifest results/run_1/checkpoints_manifest.json \\
    --output_dir         results/run_1 \\
    [--k_select_thr 0.95]

All architecture parameters (n_categories, state_dim, n_arm, latent_dim) are
read from checkpoints_manifest.json produced by 02b_train_mixvae.py.
"""

import argparse
import glob
import json
import os
import pickle
import re
import sys

import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from mpl_toolkits.axes_grid1 import make_axes_locatable

from mmidas.cpl_mixvae import cpl_mixVAE
from mmidas.eval import summarize_inference
from mmidas.utils.cluster_analysis import K_selection
from mmidas.utils.data_tools import load_data, get_loaders


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Evaluate MMIDAS checkpoints and select optimal model_order."
    )

    p.add_argument(
        "--anndata_h5ad", required=True,
        help="Path to the .h5ad file used for training (output of 01_data_prep.py)."
    )
    p.add_argument(
        "--checkpoints_manifest", required=True,
        help=(
            "Path to checkpoints_manifest.json produced by 02b_train_mixvae.py. "
            "Architecture parameters and checkpoint paths are read from this file."
        )
    )
    p.add_argument(
        "--output_dir", required=True,
        help=(
            "Directory for all outputs: summary pickle, figures, and "
            "evaluation_results.json. Created if it does not exist."
        )
    )
    p.add_argument(
        "--k_select_thr", type=float, default=0.95,
        help=(
            "Minimum average consensus threshold for K_selection(). "
            "The smallest model_order whose consensus exceeds this is selected. "
            "Default: 0.95."
        )
    )
    p.add_argument(
        "--batch_size", type=int, default=5000,
        help="Batch size for inference DataLoader. Default: 5000."
    )
    p.add_argument(
        "--seed", type=int, default=0,
        help="Random seed for train/test split (must match 02b). Default: 0."
    )

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------

def _plot_consensus_bubble(outcome, model_order, out_path):
    """Bubble-plot of the arm1-vs-arm2 consensus matrix."""
    armT_vs_armE = outcome["armA_vs_armB"][-1]
    # Use the actual matrix size — may differ from model_order if the
    # selected checkpoint's pruning state doesn't match exactly.
    k = armT_vs_armE.shape[0]
    fig, axs = plt.subplots(1, 1, figsize=(10, 10))
    mtx = armT_vs_armE / np.max(armT_vs_armE) if np.max(armT_vs_armE) > 0 else armT_vs_armE
    for l in range(k):
        for col in range(k):
            axs.add_patch(
                plt.Circle([col, l], radius=mtx[l, col], color="Navy")
            )
    axs.set_xlim([-0.5, k])
    axs.set_ylim([-0.5, k - 0.5])
    axs.invert_yaxis()
    axs.set_yticks([])
    axs.set_xticks([])
    axs.set_xlabel("arm 1", fontsize=26)
    axs.set_ylabel("arm 2", fontsize=26)
    plt.title(f"|c|= {k}", fontsize=24)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close("all")


def _plot_norm_consensus(outcome, model_order, out_path):
    """Heatmap of the normalised consensus matrix."""
    cons_mat = outcome["consensus"][0]
    fig, ax = plt.subplots(figsize=(10, 10))
    im = ax.imshow(cons_mat, cmap="rocket", vmin=0, vmax=1)
    ax.set_xlabel("T Categories (arm 1)", fontsize=30, labelpad=15)
    ax.set_ylabel("T Categories (arm 2)", fontsize=30, labelpad=15)
    ax.set_xticks([])
    ax.set_yticks([])
    plt.title(f"|c|= {cons_mat.shape[0]}", fontsize=36, pad=10)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cbar = fig.colorbar(im, cax=cax)
    for t in cbar.ax.get_yticklabels():
        t.set_fontsize(20)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close("all")
    return float(np.mean(np.diag(cons_mat)))


def _plot_state_scatter(outcome, n_arm, state_dim, model_order, output_dir):
    """Scatter plots of the continuous state variable, one per arm."""
    paths = []
    for arm in range(n_arm):
        fig = plt.figure(figsize=(5, 5))
        mu = outcome["state_mu"][arm]

        if state_dim == 1:
            ax = fig.add_subplot(1, 1, 1)
            ax.hist(mu[:, 0], bins=40)
            ax.set_xlabel("s_0")
        elif state_dim == 2:
            ax = fig.add_subplot(1, 1, 1)
            ax.scatter(mu[:, 0], mu[:, 1], s=3, alpha=0.5)
            ax.set_xlabel(r"$s_{T_1}$", fontsize=18)
            ax.set_ylabel(r"$s_{T_2}$", fontsize=18)
        else:
            ax = fig.add_subplot(1, 1, 1, projection="3d")
            ax.scatter(mu[:, 0], mu[:, 1], mu[:, 2], s=3, alpha=0.5)
            ax.set_xlabel("s_0")
            ax.set_ylabel("s_1")
            ax.set_zlabel("s_2")

        ax.set_title(f"Continuous State (arm {arm + 1})", fontsize=16, pad=12)
        fig.tight_layout()
        out_path = os.path.join(output_dir, f"state_mu_K_{model_order}_arm_{arm}.png")
        fig.savefig(out_path, dpi=150)
        plt.close("all")
        paths.append(out_path)
    return paths


# ---------------------------------------------------------------------------
# Checkpoint resolution
# ---------------------------------------------------------------------------

def resolve_checkpoint(checkpoints_manifest, n_categories, model_order):
    """Return the checkpoint path retaining `model_order` categories.

    Checkpoint naming written by 02b_train_mixvae.py:
        cpl_mixVAE_model_before_pruning_<ts>.pth        -> n_categories retained
        cpl_mixVAE_model_after_pruning_<round>_<ts>.pth -> n_categories - round

    so the pruning round we want is ``n_categories - model_order``. Round 0
    means the *before*-pruning checkpoint: there is no "after_pruning_0" file,
    and treating that miss as "no match" is what previously sent this function
    into a fallback that returned an unrelated checkpoint.

    Returns None when the model directory holds no usable checkpoint.
    """
    n_pruned_target = n_categories - model_order
    model_dir = os.path.join(os.path.dirname(checkpoints_manifest), "model")

    before_pruning = sorted(
        glob.glob(os.path.join(model_dir, "cpl_mixVAE_model_before_pruning_*.pth"))
    )
    # Map pruning round -> checkpoint path, keyed numerically. Sorting these
    # paths as strings puts round 9 after round 14, so never rely on sorted().
    after_pruning = {}
    for path in glob.glob(
        os.path.join(model_dir, "cpl_mixVAE_model_after_pruning_*.pth")
    ):
        m = re.search(r"after_pruning_(\d+)_", os.path.basename(path))
        if m:
            after_pruning[int(m.group(1))] = path

    if n_pruned_target <= 0:
        # Nothing was pruned -- this is the before-pruning checkpoint.
        if before_pruning:
            return before_pruning[-1]
    elif n_pruned_target in after_pruning:
        return after_pruning[n_pruned_target]

    if after_pruning:
        # No exact match: take the numerically closest pruning round so the
        # reported model_order stays as near the requested one as possible.
        closest = min(after_pruning, key=lambda r: (abs(r - n_pruned_target), r))
        print(
            f"  WARNING: no checkpoint for pruning round {n_pruned_target}; "
            f"using the closest available round {closest}."
        )
        return after_pruning[closest]

    if before_pruning:
        return before_pruning[-1]

    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Read architecture parameters from the training manifest
    # ------------------------------------------------------------------
    print(f"Reading training manifest: {args.checkpoints_manifest}")
    with open(args.checkpoints_manifest) as fh:
        manifest = json.load(fh)

    n_categories  = manifest["n_categories"]
    state_dim     = manifest["state_dim"]
    n_arm         = manifest["n_arm"]
    latent_dim    = manifest["latent_dim"]
    fc_dim        = manifest.get("fc_dim",        100)
    x_drop        = manifest.get("x_drop",        0.0)
    s_drop        = manifest.get("s_drop",        0.0)
    temp          = manifest.get("temp",          1.0)
    tau           = manifest.get("tau",           0.1)
    beta          = manifest.get("beta",          1.0)
    lam           = manifest.get("lam",           1.0)
    lam_pc        = manifest.get("lam_pc",        1.0)
    hard          = manifest.get("hard",          False)
    variational   = manifest.get("variational",   True)
    training_mode = manifest.get("training_mode", "MSE")
    n_gene        = manifest["n_gene"]
    checkpoints   = manifest["checkpoints"]

    if not checkpoints:
        print("ERROR: No checkpoints listed in manifest.", file=sys.stderr)
        sys.exit(1)

    print(
        f"  n_categories={n_categories}, state_dim={state_dim}, "
        f"n_arm={n_arm}, latent_dim={latent_dim}, n_gene={n_gene}"
    )
    print(f"  {len(checkpoints)} checkpoint(s) to evaluate.")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"\nLoading data: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x = data["log1p"]
    print(f"  {x.shape[0]} cells × {x.shape[1]} genes")

    # Build loaders — use the same split as training so train_loader
    # reflects the training distribution for summarize_inference.
    _, train_loader, _, _, _, _ = get_loaders(
        x=x,
        batch_size=args.batch_size,
        additional_val=False,
        seed=args.seed,
    )
    # All-data loader (no shuffling, no dropping) for final inference
    from torch.utils.data import DataLoader, TensorDataset
    import torch
    all_tensor = torch.FloatTensor(x)
    all_idx    = torch.FloatTensor(range(x.shape[0]))
    all_loader = DataLoader(
        TensorDataset(all_tensor, all_idx),
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )

    # ------------------------------------------------------------------
    # Initialise model (weights loaded per-checkpoint inside summarize_inference)
    # ------------------------------------------------------------------
    mixvae = cpl_mixVAE(saving_folder=args.output_dir, device=None)
    mixvae.init_model(
        n_categories=n_categories,
        state_dim=state_dim,
        input_dim=n_gene,
        lowD_dim=latent_dim,
        fc_dim=fc_dim,
        x_drop=x_drop,
        s_drop=s_drop,
        n_arm=n_arm,
        temp=temp,
        tau=tau,
        beta=beta,
        lam=lam,
        lam_pc=lam_pc,
        hard=hard,
        variational=variational,
        mode=training_mode,
    )
    mixvae.variational = False   # match notebook convention for evaluation

    # ------------------------------------------------------------------
    # Step 1: summarize_inference over all checkpoints → summary pickle
    # ------------------------------------------------------------------
    summary_pickle = os.path.join(
        args.output_dir,
        f"summary_performance_K_{n_categories}_narm_{n_arm}.p",
    )
    print(f"\nRunning summarize_inference over {len(checkpoints)} checkpoint(s) ...")
    summary_dict = summarize_inference(
        cpl_mixVAE=mixvae,
        files=checkpoints,
        data=train_loader,
        saving_folder=args.output_dir,
    )
    print(f"  Summary saved to: {summary_pickle}")

    # ------------------------------------------------------------------
    # Step 2: K_selection
    # ------------------------------------------------------------------
    print(f"\nRunning K_selection (thr={args.k_select_thr}) ...")
    # K_selection returns a 4-tuple:
    #   (num_pruned_sorted, recon_sorted, consensus_sorted, K)
    # K is the selected model_order (number of retained categories), or None.
    _, _, _, model_order = K_selection(
        data_dict=summary_dict,
        num_category=n_categories,
        n_arm=n_arm,
        thr=args.k_select_thr,
    )

    # Whether any checkpoint actually met k_select_thr. This is the single most
    # important fact for the human-review step, so it is recorded in
    # evaluation_results.json rather than left in the task log.
    k_selection_met_threshold = model_order is not None
    k_selection_suggested = int(model_order) if model_order is not None else None

    if model_order is None:
        print(
            f"WARNING: K_selection could not find a model meeting thr={args.k_select_thr}. "
            f"Falling back to the checkpoint with the most retained categories."
        )
        # NOTE: summary_dict["num_pruned"] holds the number of *retained*
        # categories per checkpoint (upstream misnomer -- see
        # mmidas.eval.summarize_inference, which appends len(nprune_indx)).
        model_order = int(max(summary_dict["num_pruned"]))

    model_order = int(model_order)
    print(f"  Selected model_order = {model_order}")

    # ------------------------------------------------------------------
    # Step 3: Identify selected model checkpoint
    # ------------------------------------------------------------------
    selected_model = resolve_checkpoint(
        args.checkpoints_manifest, n_categories, model_order
    )
    if selected_model is None:
        print("ERROR: Could not find a suitable model checkpoint.", file=sys.stderr)
        sys.exit(1)

    print(f"  Selected checkpoint: {selected_model}")

    # ------------------------------------------------------------------
    # Step 4: Final inference on the full dataset with the selected model
    # ------------------------------------------------------------------
    print("\nRunning final inference on full dataset ...")
    outcome = summarize_inference(
        cpl_mixVAE=mixvae,
        files=[selected_model],
        data=all_loader,
        saving_folder="",   # don't overwrite summary pickle
    )

    # Override model_order with the actual number of retained categories
    # in the loaded checkpoint — K_selection may disagree with the checkpoint
    # when the exact pruning round can't be matched.
    actual_model_order = int(len(outcome["nprune_indx"]))
    if actual_model_order != model_order:
        print(
            f"  NOTE: K_selection suggested model_order={model_order}, "
            f"but selected checkpoint has {actual_model_order} retained "
            f"categories. Using {actual_model_order}."
        )
    model_order = actual_model_order

    # ------------------------------------------------------------------
    # Step 4b: How many categories does the model actually use?
    #
    # model_order counts categories that survived pruning, not categories that
    # any cell is assigned to. A model whose discrete latent has collapsed keeps
    # a high model_order while routing every cell into a handful of categories,
    # so this is the number that says whether the run is usable.
    # ------------------------------------------------------------------
    predicted_label = outcome["pred_label"][0]        # (n_arm, n_cells), 1-indexed
    populated_per_arm = [
        int(len(np.unique(predicted_label[arm]))) for arm in range(n_arm)
    ]
    max_populated = max(populated_per_arm)
    print("\nCategory occupancy:")
    for arm, n_pop in enumerate(populated_per_arm):
        print(f"  arm {arm}: {n_pop} of {model_order} retained categories are non-empty")

    collapse_warning = None
    if max_populated < 0.5 * model_order:
        collapse_warning = (
            f"Only {max_populated} of {model_order} retained categories are "
            f"populated (arms: {populated_per_arm}). The discrete latent has "
            f"likely collapsed -- model_order is not a usable estimate of the "
            f"number of cell types and downstream Analyze results will be "
            f"dominated by empty categories."
        )
        print(f"  WARNING: {collapse_warning}")

    # ------------------------------------------------------------------
    # Step 5: Figures
    # ------------------------------------------------------------------
    print("\nSaving figures ...")
    figures = []

    bubble_path = os.path.join(args.output_dir, f"consensus_T1_vs_T2_K_{model_order}.png")
    _plot_consensus_bubble(outcome, model_order, bubble_path)
    figures.append(bubble_path)
    print(f"  {bubble_path}")

    norm_path = os.path.join(args.output_dir, f"norm_consensus_T1_vs_T2_K_{model_order}.png")
    avg_consensus = _plot_norm_consensus(outcome, model_order, norm_path)
    figures.append(norm_path)
    print(f"  {norm_path}  (avg_consensus={avg_consensus:.4f})")

    state_paths = _plot_state_scatter(outcome, n_arm, state_dim, model_order, args.output_dir)
    figures.extend(state_paths)
    for p in state_paths:
        print(f"  {p}")

    # ------------------------------------------------------------------
    # Step 6: evaluation_results.json
    # ------------------------------------------------------------------
    results = {
        "model_order":      model_order,
        "n_categories":     n_categories,
        "n_arm":            n_arm,
        "state_dim":        state_dim,
        "latent_dim":       latent_dim,
        "n_gene":           n_gene,
        "selected_model":   selected_model,
        "summary_pickle":   summary_pickle,
        "avg_consensus":    avg_consensus,
        "k_select_thr":     args.k_select_thr,
        # Review fields: a run can produce a plausible-looking model_order while
        # failing both of these, so they travel with the hand-off file.
        "k_selection_met_threshold":        k_selection_met_threshold,
        "k_selection_suggested_model_order": k_selection_suggested,
        "n_populated_categories":          max_populated,
        "n_populated_categories_per_arm":  populated_per_arm,
        "collapse_warning":                collapse_warning,
        "figures":          figures,
    }
    results_path = os.path.join(args.output_dir, "evaluation_results.json")
    with open(results_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print(f"\nEvaluation complete.")
    print(f"  model_order            = {model_order}")
    print(f"  n_populated_categories = {max_populated}")
    print(f"  avg_consensus          = {avg_consensus:.4f}")
    print(f"  k_select_thr met       = {k_selection_met_threshold}")
    print(f"  Results JSON   : {results_path}")

    if not k_selection_met_threshold or collapse_warning:
        print(
            "\nREVIEW REQUIRED: this model did not pass automated selection. "
            "Do not launch MMIDAS_Analyze on it without addressing the warnings "
            "above."
        )


if __name__ == "__main__":
    sys.exit(main())
