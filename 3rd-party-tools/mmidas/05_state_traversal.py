#!/usr/bin/env python3
"""
05_state_traversal.py — MMIDAS State Traversal Figures

Uses the pickles written by 03c_traversal_prep.py to produce:

  A. State-space scatter (always):
       state/s_mu/state_mu_arm_{arm}_c_{cc}.png  — per selected category
       All cells plotted in gray, target category in color, traversal path
       overlaid in black.

  B. Per-pathway box plots (only when g_subset is non-empty):
       state/signaling_pathways/s_pc_path_{ig+1}_K_{K}.png
       For each KEGG pathway: box plots of gene variation V_s(g) across
       selected categories with scatter overlay.

  C. All-pathway sorted box plots per selected category (only when g_subset
     is non-empty):
       state/pathway_summary/pathway_summary_c_{cc}.png

Usage:
  python3 05_state_traversal.py \\
    --anndata_h5ad            Mouse_ALM-VISp_cpm.h5ad \\
    --checkpoints_manifest    results/run_1/checkpoints_manifest.json \\
    --evaluation_results_json results/run_1/evaluation_results.json \\
    --traversal_manifest      results/run_1/traversal_manifest.json \\
    [--arm              0]               # arm to use for s_mu scatter
    [--n_selected_cats  10]              # limit categories; 0 = all active
    [--batch_size       5000]
    [--seed             0]
"""

import argparse
import json
import os
import pickle
import sys
import warnings

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import torch
from torch.utils.data import DataLoader, TensorDataset

from mmidas.cpl_mixvae import cpl_mixVAE
from mmidas.eval import summarize_inference
from mmidas.utils.data_tools import load_data

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate state traversal figures from 03c traversal pickles."
    )
    p.add_argument("--anndata_h5ad", required=True,
                   help="Path to .h5ad data file (output of 01_data_prep.py).")
    p.add_argument("--checkpoints_manifest", required=True,
                   help="Path to checkpoints_manifest.json from 02b_train_mixvae.py.")
    p.add_argument("--evaluation_results_json", required=True,
                   help="Path to evaluation_results.json from 03a_evaluate.py.")
    p.add_argument("--traversal_manifest", required=True,
                   help="Path to traversal_manifest.json from 03c_traversal_prep.py.")
    p.add_argument("--arm", type=int, default=0,
                   help="Which arm to use for state-space scatter plots. Default: 0.")
    p.add_argument("--n_selected_cats", type=int, default=10,
                   help=(
                       "Number of active categories to generate pathway/scatter plots "
                       "for (sorted by category index). 0 = all active. Default: 10."
                   ))
    p.add_argument("--batch_size", type=int, default=5000,
                   help="Batch size for model inference. Default: 5000.")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed. Default: 0.")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _label_for_cat(cc: int, predicted_label_arm: np.ndarray,
                   cluster_labels: np.ndarray) -> str:
    """Return a short annotation string for category cc (1-indexed)."""
    from statistics import mode as stat_mode
    cells = np.where(predicted_label_arm == cc)[0]
    if len(cells) == 0:
        return f"Cat_{cc}"
    dom_ttype = stat_mode(cluster_labels[cells])
    # Truncating at the first space collapses distinct t-types to the same
    # word -- "L6 IT Car3" and "L6 CT Nxph2" both became "L6", so a ten-panel
    # figure could be labelled L6, Pvalb, L6, L6, L4, Vip, L6, L6, Pvalb, Pvalb
    # with no way to tell the categories apart. Keep a short form of the
    # dominant t-type but always qualify it with the category number.
    space = dom_ttype.find(" ")
    short = dom_ttype[:space] if space > 0 else dom_ttype[:12]
    return f"c{cc} {short}"


def _default_palette(n: int):
    """Return n distinct hex colors that stand out against the grey background.

    Not tab10: its 8th entry is grey (#7f7f7f), which is invisible against the
    #dbdbde background used for unselected cells -- the 8th category's panel
    looked monochrome and empty. husl is evenly spaced around the hue circle
    and never returns a neutral.
    """
    return sns.color_palette("husl", n_colors=n).as_hex()


# ---------------------------------------------------------------------------
# Plot A: state-space scatter
# ---------------------------------------------------------------------------

def _plot_state_scatter(
    state_mu_all: np.ndarray,  # (n_cells, state_dim) for one arm
    s_travers: list,           # list of (n_arm, n_steps, state_dim) per category
    c_cat: np.ndarray,         # 1-indexed active categories
    selected_c: list,          # 1-indexed subset to plot
    predicted_label_arm: np.ndarray,  # (n_cells,) 1-indexed
    cluster_labels: np.ndarray,
    colors: list,
    out_dir: str,
    arm: int,
    model_order: int,
):
    s_mu_dir = os.path.join(out_dir, "state", "s_mu")
    os.makedirs(s_mu_dir, exist_ok=True)

    for i_c, cc in enumerate(selected_c):
        ic = np.where(c_cat == cc)[0][0]
        idx_c = np.where(predicted_label_arm == cc)[0]

        fig = plt.figure(figsize=(5, 5))
        sns.set_theme()
        sns.set(rc={"axes.facecolor": "#ffffff"})
        ax = fig.add_subplot(1, 1, 1)

        # All cells in gray
        ax.scatter(state_mu_all[:, 0], state_mu_all[:, 1],
                   color="#dbdbde", s=2, alpha=0.7)
        # Target category highlighted
        if len(idx_c) > 0:
            ax.scatter(state_mu_all[idx_c, 0], state_mu_all[idx_c, 1],
                       color=colors[i_c], s=6, alpha=1.0)
        # Traversal path
        trav = s_travers[ic][arm]  # (n_steps, state_dim)
        ax.scatter(trav[:, 0], trav[:, 1], color="black", s=0.5, alpha=1.0)

        label = _label_for_cat(cc, predicted_label_arm, cluster_labels)
        ax.set_xlabel(r"$s_{T_1}$", fontsize=16)
        ax.set_ylabel(r"$s_{T_2}$", fontsize=16)
        ax.set_title(label, fontsize=16, pad=8)
        ax.xaxis.set_tick_params(labelsize=9)
        ax.yaxis.set_tick_params(labelsize=9)
        ax.grid(False)
        for sp in ax.spines.values():
            sp.set_color("#2d2d2d")
        fig.tight_layout()
        out_path = os.path.join(s_mu_dir, f"state_mu_arm_{arm}_c_{cc}.png")
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close("all")

    print(f"  Saved {len(selected_c)} state-space scatter(s) → {s_mu_dir}")


# ---------------------------------------------------------------------------
# Plot B: per-pathway box plots
# ---------------------------------------------------------------------------

def _plot_pathway_boxes(
    g_var_mean: np.ndarray,  # (n_cats, n_arm, 1, n_gene)
    g_subset: list,
    pathways: list,
    c_cat: np.ndarray,
    selected_c: list,
    c_labels: list,           # annotation strings for selected_c
    colors: list,
    out_dir: str,
    arm: int,
    model_order: int,
) -> list:
    sig_dir = os.path.join(out_dir, "state", "signaling_pathways")
    os.makedirs(sig_dir, exist_ok=True)
    produced = []

    boxprops    = dict(color="black", linewidth=1)
    medianprops = dict(color="black", linewidth=1.5)

    for ig, g_s in enumerate(g_subset):
        sig_path = pathways[ig]
        g_s = np.array(g_s)

        vals, names, xs = [], [], []
        med_value = []
        rng = np.random.default_rng(0)

        for i_c, cc in enumerate(selected_c):
            c_idx = np.where(c_cat == cc)[0][0]
            g_travers = g_var_mean[c_idx, arm, 0, g_s]
            med_value.append(np.median(g_travers))
            vals.append(g_travers)
            names.append(c_labels[i_c])
            xs.append(rng.normal(i_c + 1, 0.04, g_travers.shape[0]))

        plt.close("all")
        plt.figure(figsize=[max(8, len(selected_c) * 1.2), 5])
        overall_med = np.median(med_value)
        plt.plot(
            np.linspace(0.5, len(selected_c) + 0.5, 10),
            overall_med * np.ones(10),
            linestyle="--", linewidth=0.8, color="gray",
        )
        plt.boxplot(vals, labels=names, showfliers=False,
                    boxprops=boxprops, medianprops=medianprops)
        for x, val, c in zip(xs, vals, colors):
            plt.scatter(x, val, alpha=0.4, color=c)
        plt.title(sig_path, fontsize=22, pad=16)
        plt.ylabel(r"$V_{s}(g)$", fontsize=20, labelpad=8)
        plt.xticks(fontsize=14, rotation=90)
        plt.ylim([0, 1.5])
        plt.yticks([0, 0.5, 1, 1.5], fontsize=14)
        plt.tight_layout()
        out_path = os.path.join(sig_dir, f"s_pc_path_{ig + 1}_K_{model_order}.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        produced.append(out_path)

    print(f"  Saved {len(produced)} pathway box plot(s) → {sig_dir}")
    return produced


# ---------------------------------------------------------------------------
# Plot C: all-pathway sorted box plots per category
# ---------------------------------------------------------------------------

def _plot_pathway_summary(
    g_var_mean: np.ndarray,
    g_subset: list,
    pathways: list,
    c_cat: np.ndarray,
    selected_c: list,
    c_labels: list,
    colors: list,
    out_dir: str,
    arm: int,
    model_order: int,
) -> list:
    summary_dir = os.path.join(out_dir, "state", "pathway_summary")
    os.makedirs(summary_dir, exist_ok=True)
    produced = []

    boxprops    = dict(color="black", linewidth=1)
    medianprops = dict(color="black", linewidth=1.5)

    for i_c, cc in enumerate(selected_c):
        c_idx = np.where(c_cat == cc)[0][0]
        vals, names, path_med = [], [], []

        for ig, g_s in enumerate(g_subset):
            g_s = np.array(g_s)
            g_travers = g_var_mean[c_idx, arm, 0, g_s]
            path_med.append(np.median(g_travers))
            vals.append(g_travers)
            names.append(pathways[ig])

        path_med = np.array(path_med)
        sort_idx = np.argsort(path_med)
        vals  = [vals[ii]  for ii in sort_idx]
        names = [names[ii] for ii in sort_idx]

        plt.close("all")
        plt.figure(figsize=[max(10, len(g_subset) * 0.7), 4])
        plt.boxplot(vals, labels=names, showfliers=False,
                    boxprops=boxprops, medianprops=medianprops)
        for iv, val in enumerate(vals):
            rng = np.random.default_rng(iv)
            jitter = rng.normal(0, 0.04, len(val))
            plt.scatter((iv + 1) * np.ones(len(val)) + jitter, val,
                        alpha=0.4, color=colors[i_c])
        plt.title(c_labels[i_c], fontsize=18)
        plt.ylabel(r"$V_{s}(g)$", fontsize=16)
        plt.xticks(fontsize=11, rotation=90)
        plt.ylim([0, 1.5])
        plt.tight_layout()
        out_path = os.path.join(summary_dir, f"pathway_summary_c_{cc}.png")
        plt.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close("all")
        produced.append(out_path)

    print(f"  Saved {len(produced)} pathway summary plot(s) → {summary_dir}")
    return produced


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    np.random.seed(args.seed)

    # ------------------------------------------------------------------
    # Load manifests
    # ------------------------------------------------------------------
    with open(args.checkpoints_manifest) as fh:
        ckpt_mf = json.load(fh)
    with open(args.evaluation_results_json) as fh:
        eval_results = json.load(fh)
    with open(args.traversal_manifest) as fh:
        trav_mf = json.load(fh)

    n_categories  = ckpt_mf["n_categories"]
    state_dim     = ckpt_mf["state_dim"]
    n_arm         = ckpt_mf["n_arm"]
    latent_dim    = ckpt_mf["latent_dim"]
    fc_dim        = ckpt_mf.get("fc_dim", 100)
    x_drop        = ckpt_mf.get("x_drop", 0.0)
    s_drop        = ckpt_mf.get("s_drop", 0.0)
    temp          = ckpt_mf.get("temp", 1.0)
    tau           = ckpt_mf.get("tau", 0.1)
    beta          = ckpt_mf.get("beta", 1.0)
    lam           = ckpt_mf.get("lam", 1.0)
    lam_pc        = ckpt_mf.get("lam_pc", 1.0)
    hard          = ckpt_mf.get("hard", False)
    variational   = ckpt_mf.get("variational", True)
    training_mode = ckpt_mf.get("training_mode", "MSE")
    n_gene        = ckpt_mf["n_gene"]

    model_order    = eval_results["model_order"]
    selected_model = eval_results["selected_model"]
    out_dir        = trav_mf["output_dir"]

    arm = min(args.arm, n_arm - 1)

    print(f"model_order={model_order}, n_arm={n_arm}, arm={arm}, state_dim={state_dim}")

    # ------------------------------------------------------------------
    # Load traversal pickles
    # ------------------------------------------------------------------
    print("\nLoading traversal pickles ...")
    with open(trav_mf["traversal_pickle"], "rb") as fh:
        trav_data = pickle.load(fh)
    with open(trav_mf["state_mu_pickle"], "rb") as fh:
        smu_data = pickle.load(fh)

    g_var_mean = trav_data["V_g_mean"]   # (n_cats, n_arm, 1, n_gene)
    g_subset   = trav_data["g_subset"]   # list of gene-index arrays
    pathways   = trav_data["pathways"]   # list of pathway names
    c_cat      = trav_data["c_cat"]      # 1-indexed active categories
    s_travers  = smu_data["s_travers"]   # list of (n_arm, n_steps, state_dim)

    n_cats = len(c_cat)
    print(f"  {n_cats} active categories, {len(pathways)} KEGG pathways loaded.")

    # Categories are selected once per-cell assignments are available, further
    # down -- picking them here by index alone can select categories that no
    # cell was assigned to, which yields identical, empty figures.
    n_sel = args.n_selected_cats if args.n_selected_cats > 0 else n_cats
    n_sel = min(n_sel, n_cats)

    # ------------------------------------------------------------------
    # Data + model inference (for full per-cell state_mu)
    # ------------------------------------------------------------------
    print(f"\nLoading data: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x              = data["log1p"]
    cluster_labels = np.array([c.strip() for c in data["cluster"]])
    n_cells        = x.shape[0]

    print("Running model inference ...")
    mixvae = cpl_mixVAE(saving_folder=out_dir, device=None)
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
    mixvae.variational = False

    loader = DataLoader(
        TensorDataset(torch.FloatTensor(x), torch.FloatTensor(range(n_cells))),
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )

    outcome = summarize_inference(
        cpl_mixVAE=mixvae,
        files=[selected_model],
        data=loader,
        saving_folder="",
    )

    # state_mu shape: (n_arm, n_cells, state_dim)
    # pred_label[0] shape: (n_arm, n_cells) — 1-indexed
    state_mu_all    = outcome["state_mu"][arm]    # (n_cells, state_dim)
    predicted_label = outcome["pred_label"][0][arm]  # (n_cells,) 1-indexed

    # ------------------------------------------------------------------
    # Select subset of categories, largest first.
    #
    # "Active" (survived pruning) is not the same as "populated". Selecting the
    # first n_selected_cats by category index can land entirely on categories
    # with zero assigned cells, in which case every figure is the same plot of
    # the background scatter with a degenerate traversal path. Rank by assigned
    # cell count so the figures show the categories the model actually uses.
    # ------------------------------------------------------------------
    cell_counts = {
        int(cc): int(np.sum(predicted_label == cc)) for cc in c_cat
    }
    populated = [cc for cc in c_cat if cell_counts[int(cc)] > 0]
    print(
        f"\nCategory occupancy (arm {arm}): {len(populated)} of {n_cats} active "
        f"categories are non-empty."
    )
    if len(populated) < n_cats:
        print(
            f"  WARNING: {n_cats - len(populated)} active categories have no "
            f"cells assigned. They are excluded from the figures below."
        )

    if not populated:
        print(
            "ERROR: no active category has any cells assigned to it -- the "
            "model's discrete latent has collapsed and there is nothing to "
            "plot. Re-check the Train stage before running Analyze.",
            file=sys.stderr,
        )
        return 1

    ranked = sorted(populated, key=lambda cc: (-cell_counts[int(cc)], int(cc)))
    selected_c = [int(cc) for cc in ranked[:n_sel]]
    n_sel = len(selected_c)
    print(
        f"  Selected {n_sel} category/categories by cell count: "
        + ", ".join(f"{cc} (n={cell_counts[cc]})" for cc in selected_c)
    )

    # ------------------------------------------------------------------
    # Annotations and colors
    # ------------------------------------------------------------------
    colors = _default_palette(n_sel)
    c_labels = [
        _label_for_cat(int(cc), predicted_label, cluster_labels)
        for cc in selected_c
    ]

    # ------------------------------------------------------------------
    # Plot A: state-space scatter
    # ------------------------------------------------------------------
    print(f"\nPlot A: state-space scatter ({n_sel} categories, arm {arm}) ...")
    _plot_state_scatter(
        state_mu_all=state_mu_all,
        s_travers=s_travers,
        c_cat=c_cat,
        selected_c=selected_c,
        predicted_label_arm=predicted_label,
        cluster_labels=cluster_labels,
        colors=colors,
        out_dir=out_dir,
        arm=arm,
        model_order=model_order,
    )

    pathway_figs, summary_figs = [], []

    if g_subset:
        # ------------------------------------------------------------------
        # Plot B: per-pathway box plots
        # ------------------------------------------------------------------
        print(f"\nPlot B: per-pathway box plots ({len(pathways)} pathways) ...")
        pathway_figs = _plot_pathway_boxes(
            g_var_mean=g_var_mean,
            g_subset=g_subset,
            pathways=pathways,
            c_cat=c_cat,
            selected_c=selected_c,
            c_labels=c_labels,
            colors=colors,
            out_dir=out_dir,
            arm=arm,
            model_order=model_order,
        )

        # ------------------------------------------------------------------
        # Plot C: all-pathway sorted per-category
        # ------------------------------------------------------------------
        print(f"\nPlot C: all-pathway summary per category ({n_sel} categories) ...")
        summary_figs = _plot_pathway_summary(
            g_var_mean=g_var_mean,
            g_subset=g_subset,
            pathways=pathways,
            c_cat=c_cat,
            selected_c=selected_c,
            c_labels=c_labels,
            colors=colors,
            out_dir=out_dir,
            arm=arm,
            model_order=model_order,
        )
    else:
        print("\nNo KEGG pathways in traversal pickle — skipping pathway figures.")

    # ------------------------------------------------------------------
    # Manifest
    # ------------------------------------------------------------------
    s_mu_dir = os.path.join(out_dir, "state", "s_mu")
    scatter_figs = [
        os.path.join(s_mu_dir, f"state_mu_arm_{arm}_c_{cc}.png")
        for cc in selected_c
    ]

    manifest = {
        "output_dir":      out_dir,
        "model_order":     model_order,
        "arm":             arm,
        "n_selected_cats": n_sel,
        "selected_c":      [int(c) for c in selected_c],
        # Cell count per selected category, and how many of the model's active
        # categories are populated at all -- lets a validation step tell a real
        # traversal from figures drawn over empty categories.
        "selected_c_n_cells":       [cell_counts[int(c)] for c in selected_c],
        "n_active_categories":      int(n_cats),
        "n_populated_categories":   int(len(populated)),
        "n_pathways":      len(pathways),
        "scatter_figs":    scatter_figs,
        "pathway_figs":    pathway_figs,
        "summary_figs":    summary_figs,
    }
    mani_path = os.path.join(out_dir, "state_traversal_manifest.json")
    with open(mani_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    total = len(scatter_figs) + len(pathway_figs) + len(summary_figs)
    print(f"\nState traversal complete — {total} total figures.")
    print(f"  Manifest: {mani_path}")


if __name__ == "__main__":
    sys.exit(main())
