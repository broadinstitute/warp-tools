#!/usr/bin/env python3
"""
04_clusterability.py — MMIDAS Clusterability Figures

Loads the classification pickles written by 03b_classify.py and produces:

  1. classAcc_RF_K_{K}.png
       Grouped bar chart: RF classification accuracy for three label sets
       (t-types, ConsType arm 0, ConsType arm 1) × two feature embeddings
       (Linear/PCA, Non-Linear/MMIDAS lowD). Error bars = 1 SD across k-folds.

  2. SC_K_{K}_{date}.png
       Per-category silhouette scores (sorted). One curve per arm for ConsType
       categories; one reference curve for t-types (PCA features).

  3. Confusion-matrix figures (Blues heatmaps, one per combination):
       conf_Ttype_pc.png                   — t-type × PCA
       conf_Ttype_lowD_arm_{arm}.png        — t-type × MMIDAS arm
       conf_ConsType_pc_arm_{arm}.png       — ConsType × PCA arm
       conf_ConsType_lowD_arm_{arm}.png     — ConsType × MMIDAS arm

Usage:
  python3 04_clusterability.py \\
    --classify_manifest  results/run_1/classify_manifest.json \\
    [--output_dir        results/run_1]        # defaults to same dir as manifest
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
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load(path: str) -> dict:
    with open(path, "rb") as fh:
        return pickle.load(fh)


def _conf_heatmap(mat: np.ndarray, ylabel: str, out_path: str) -> None:
    """Save a Blues confusion-matrix heatmap (zeros masked).

    `mat` arrives as raw integer counts from the classification pickles, but the
    colour scale below is fixed to [0, 1]. Plotting counts directly saturates
    every non-zero cell to the darkest blue, so a single misassigned cell looks
    exactly as strong as a correct diagonal entry of 860 -- the figure collapses
    to black-and-white noise and carries no information.
    Row-normalise first, matching the reference analysis
    (notebooks/3_evaluation.ipynb divides each row by its sum before plotting),
    so each cell reads as "fraction of this true class predicted as that".
    """
    eps = 1e-3
    mat = mat.copy().astype(float)
    row_sum = mat.sum(axis=1, keepdims=True)
    mat = np.divide(mat, row_sum, out=np.zeros_like(mat), where=row_sum != 0)
    mat[mat < eps] = 0.0
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.set_theme()
    sns.set(rc={"axes.facecolor": "white"})
    sns.heatmap(
        mat,
        xticklabels=[],
        yticklabels=[],
        vmin=0,
        vmax=1,
        ax=ax,
        cmap="Blues",
        cbar_kws={"shrink": 1},
        mask=(mat == 0),
    )
    ax.set_ylabel(ylabel, fontsize=22, labelpad=12)
    ax.set_xlabel("Prediction", fontsize=22, labelpad=12)
    for spine_fn in (ax.axhline, ax.axvline):
        pass  # border lines below
    ax.axhline(y=0,          color="#2d2d2d", linewidth=1)
    ax.axhline(y=mat.shape[0], color="#2d2d2d", linewidth=1)
    ax.axvline(x=0,          color="#2d2d2d", linewidth=1)
    ax.axvline(x=mat.shape[1], color="#2d2d2d", linewidth=1)
    fig.tight_layout()
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Generate clusterability figures from 03b classification pickles."
    )
    p.add_argument("--classify_manifest", required=True,
                   help="Path to classify_manifest.json from 03b_classify.py.")
    p.add_argument("--output_dir", default="",
                   help=(
                       "Directory for output figures. "
                       "Defaults to the directory containing classify_manifest.json."
                   ))
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    with open(args.classify_manifest) as fh:
        mf = json.load(fh)

    clustering_dir = mf["clustering_dir"]
    model_order    = mf["model_order"]
    n_ttype        = mf["n_ttype"]
    n_pca          = mf["n_pca"]
    latent_dim     = mf["latent_dim"]
    n_arm          = mf["n_arm"]
    k_fold         = mf["k_fold"]
    date_tag       = mf["date_tag"]

    out_dir = args.output_dir or os.path.dirname(args.classify_manifest)
    os.makedirs(out_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load pickles
    # ------------------------------------------------------------------
    print("Loading classification pickles ...")

    def _p(name):
        return os.path.join(clustering_dir, name)

    # 1. Ttype × PCA  (global, no arm)
    ttype_pc = _load(_p(
        f"Ttype_classification_K_{n_ttype}_nFeature_{n_pca}_{date_tag}.p"
    ))
    acc_T_pc   = ttype_pc["acc_T_adj"]   # shape (k_fold,)
    sc_T_pc    = ttype_pc["sc_T"]        # shape (1, n_ttype)
    conf_T_pc  = ttype_pc["conf_mat"]

    # Per-arm data
    acc_T_lowD,      sc_T_lowD,      conf_T_lowD      = [], [], []
    acc_cons_pc,     sc_cons_pc,     conf_cons_pc      = [], [], []
    acc_cons_lowD,   sc_cons_lowD,   conf_cons_lowD    = [], [], []

    for arm in range(n_arm):
        # Ttype × MMIDAS arm
        d = _load(_p(
            f"Ttype_classification_K_{n_ttype}_nFeature_{latent_dim}_arm_{arm}_{date_tag}.p"
        ))
        acc_T_lowD.append(d["acc_T_adj"])
        sc_T_lowD.append(d["sc_T"])
        conf_T_lowD.append(d["conf_mat"])

        # ConsType × PCA arm
        d = _load(_p(
            f"ConsType_classification_K_{model_order}_nFeature_{n_pca}_arm_{arm}_{date_tag}.p"
        ))
        acc_cons_pc.append(d["acc_T_adj"])
        sc_cons_pc.append(d["sc_T"])
        conf_cons_pc.append(d["conf_mat"])

        # ConsType × MMIDAS arm
        d = _load(_p(
            f"ConsType_classification_K_{model_order}_nFeature_{latent_dim}_arm_{arm}_{date_tag}.p"
        ))
        acc_cons_lowD.append(d["acc_T_adj"])
        sc_cons_lowD.append(d["sc_T"])
        conf_cons_lowD.append(d["conf_mat"])

    # ------------------------------------------------------------------
    # Figure 1: Bar chart — classification accuracy
    # ------------------------------------------------------------------
    print("\nFigure 1: classification accuracy bar chart ...")

    # Average lowD accuracy across arms for the t-type label set
    acc_T_lowDs = np.mean(np.vstack(acc_T_lowD), axis=0)  # (k_fold,)

    n_cond = 2  # PCA vs lowD
    label_col, embed_col, acc_col = [], [], []

    for arm in range(n_arm):
        # Newline rather than a long single line: rotated long labels overlap
        # each other once there is more than one arm.
        group_label = f"T Categories\n(arm {arm + 1})"
        label_col.extend([group_label] * k_fold)
        embed_col.extend(["Linear (PCA)"] * k_fold)
        acc_col.extend(acc_cons_pc[arm] * 100)

        label_col.extend([group_label] * k_fold)
        embed_col.extend(["Non-Linear (MMIDAS)"] * k_fold)
        acc_col.extend(acc_cons_lowD[arm] * 100)

    # Prepend t-type group so it's leftmost
    all_labels  = ["t-types"] * (n_cond * k_fold) + label_col
    all_embeds  = (["Linear (PCA)"] * k_fold + ["Non-Linear (MMIDAS)"] * k_fold) + embed_col
    all_acc     = list(acc_T_pc * 100) + list(acc_T_lowDs * 100) + acc_col

    df = pd.DataFrame({
        "Accuracy":        all_acc,
        "Cell Types":      all_labels,
        "Class_Embedding": all_embeds,
    })

    plt.close("all")
    # Wide enough for one "T Categories (arm N)" label per arm plus the
    # t-types group. At the previous width these ran together as
    # "t-typesT CategoriesT Categories".
    fig, ax = plt.subplots(figsize=(max(8, 4 + 2.5 * n_arm), 6))
    sns.reset_defaults()
    pal = sns.color_palette("Reds", n_colors=6)
    bar_colors = [pal.as_hex()[1], pal.as_hex()[4]]

    sns.barplot(
        data=df,
        x="Cell Types",
        y="Accuracy",
        hue="Class_Embedding",
        ax=ax,
        width=0.8,
        errorbar="sd",
        capsize=0.1,
        errwidth=2,
        edgecolor="black",
        palette=bar_colors,
    )

    # Add hatching on "Non-Linear" bars for accessibility
    hatch_map = {"Linear (PCA)": "", "Non-Linear (MMIDAS)": "//"}
    legend_handles = ax.get_legend().legend_handles if ax.get_legend() else []
    for container, handle in zip(ax.containers, legend_handles):
        hatch = hatch_map.get(handle.get_label(), "")
        for patch in container:
            patch.set_hatch(hatch)
        handle.set_hatch(hatch + hatch)

    ax.legend(title="Embedding", frameon=False, loc="upper left",
              fontsize=11, bbox_to_anchor=(1, 1))
    ax.set_ylim([0, 101])
    ax.set_ylabel("Accuracy (%)", fontsize=18, labelpad=10)
    ax.set_xlabel("")
    ax.xaxis.set_tick_params(labelsize=13)
    ax.yaxis.set_tick_params(labelsize=13)
    ax.xaxis.set_ticks_position("none")
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(0)
        lbl.set_horizontalalignment("center")
    ax.grid(False)
    fig.tight_layout()
    acc_png = os.path.join(out_dir, f"classAcc_RF_K_{model_order}.png")
    fig.savefig(acc_png, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {acc_png}")

    # ------------------------------------------------------------------
    # Figure 2: Silhouette scores
    # ------------------------------------------------------------------
    print("Figure 2: silhouette scores ...")

    pal_r = sns.color_palette("Reds", n_colors=6)
    mks = ["^", "v", ">", "<"]

    plt.close("all")
    plt.figure(figsize=(5, 6), dpi=100)

    # The classification pickles store sc_T with a leading singleton axis, i.e.
    # shape (1, n_category). Sorting that 2-D array and taking len() gives 1, so
    # every point ends up plotted at x=1 as a separate broadcast line -- each
    # inheriting the same label, which is what produced a figure thousands of
    # pixels tall with one legend entry per category. Flatten first.
    # t-type reference line (PCA)
    sc_T_sorted = np.sort(np.ravel(sc_T_pc))
    n_sc = len(sc_T_sorted)
    plt.plot(np.arange(n_sc) + 1, sc_T_sorted, linewidth=2, linestyle="--",
             label="t-types (PCA)", color=pal_r[1])

    # ConsType per arm
    for arm in range(n_arm):
        sc_arm = np.sort(np.ravel(sc_cons_lowD[arm]))
        n_sc_arm = len(sc_arm)
        plt.plot(
            np.arange(n_sc_arm) + 1,
            sc_arm,
            linewidth=1.5,
            marker=mks[arm % len(mks)],
            markersize=4,
            label=f"T Category (arm {arm + 1})",
            color=pal_r[2 + arm],
        )

    plt.xlabel("# Cell Type", fontsize=18, labelpad=12)
    plt.xticks(fontsize=14)
    plt.ylabel("Silhouette Score", fontsize=18, labelpad=10)
    plt.yticks(fontsize=14)
    plt.ylim([-0.25, 1.0])
    plt.legend(fontsize=11, loc="lower right", frameon=False)
    plt.tight_layout()
    sc_png = os.path.join(out_dir, f"SC_K_{model_order}_{date_tag}.png")
    plt.savefig(sc_png, dpi=300, bbox_inches="tight")
    plt.close("all")
    print(f"  Saved: {sc_png}")

    # ------------------------------------------------------------------
    # Figure 3: Confusion matrices
    # ------------------------------------------------------------------
    print("Figure 3: confusion matrices ...")

    _conf_heatmap(
        conf_T_pc,
        "T-type",
        os.path.join(out_dir, "conf_Ttype_pc.png"),
    )

    for arm in range(n_arm):
        _conf_heatmap(
            conf_T_lowD[arm],
            "T-type",
            os.path.join(out_dir, f"conf_Ttype_lowD_arm_{arm}.png"),
        )
        _conf_heatmap(
            conf_cons_pc[arm],
            "T Category",
            os.path.join(out_dir, f"conf_ConsType_pc_arm_{arm}.png"),
        )
        _conf_heatmap(
            conf_cons_lowD[arm],
            "T Category",
            os.path.join(out_dir, f"conf_ConsType_lowD_arm_{arm}.png"),
        )

    # ------------------------------------------------------------------
    # Manifest
    # ------------------------------------------------------------------
    produced = [acc_png, sc_png,
                os.path.join(out_dir, "conf_Ttype_pc.png")]
    for arm in range(n_arm):
        produced += [
            os.path.join(out_dir, f"conf_Ttype_lowD_arm_{arm}.png"),
            os.path.join(out_dir, f"conf_ConsType_pc_arm_{arm}.png"),
            os.path.join(out_dir, f"conf_ConsType_lowD_arm_{arm}.png"),
        ]

    manifest = {
        "output_dir":  os.path.abspath(out_dir),
        "model_order": model_order,
        "n_ttype":     n_ttype,
        "n_arm":       n_arm,
        "figures":     produced,
    }
    mani_path = os.path.join(out_dir, "clusterability_manifest.json")
    with open(mani_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nClusterability figures complete — {len(produced)} figures.")
    print(f"  Manifest: {mani_path}")


if __name__ == "__main__":
    sys.exit(main())
