#!/usr/bin/env python3
"""
03c_traversal_prep.py — MMIDAS State Traversal Data Preparation

Pre-computes all data needed by 05_state_traversal.py:

  1. state/traversal_pca_K_{model_order}.pickle
       keys: V_g_mean, V_g_std, g_subset, c_cat, pathways
       V_g_mean[cat_idx, arm, 0, gene] = per-gene std of reconstructed expression
       along a PCA traversal of the per-category state distribution.

  2. state/state_mu_pca_K_{model_order}.pickle
       keys: s_mu, s_travers, c_cat
       s_mu[cat_idx]    = mean state for each category (n_arm, state_dim)
       s_travers[cat_idx] = traversal path in state space (n_arm, n_steps, state_dim)

  3. taxonomy_order_K_{model_order}.npy
       Array of category indices (0-based) ordered by their dominant t-type's
       position in the taxonomy tree. Used to assign T_N labels in notebook 5.
       Falls back to alphabetical t-type order if --htree_file is not supplied.

Usage:
  python3 03c_traversal_prep.py \\
    --anndata_h5ad            Mouse_ALM-VISp_cpm.h5ad \\
    --checkpoints_manifest    results/run_1/checkpoints_manifest.json \\
    --evaluation_results_json results/run_1/evaluation_results.json \\
    --output_dir              results/run_1 \\
    [--kegg_toml              KEGG/KEGG.toml] \\
    [--htree_file             data/tree_Mouse_ALM-VISp_2018.csv] \\
    [--n_traversal_steps      50]

For a fast smoke test use: --n_traversal_steps 5
"""

import argparse
import json
import os
import pickle
import sys

import matplotlib
matplotlib.use("Agg")

import numpy as np
import torch
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader, TensorDataset

from mmidas.cpl_mixvae import cpl_mixVAE
from mmidas.eval import summarize_inference
from mmidas.utils.data_tools import load_data


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Pre-compute state traversal data for 05_state_traversal.py."
    )
    p.add_argument("--anndata_h5ad", required=True,
                   help="Path to the .h5ad file (output of 01_data_prep.py).")
    p.add_argument("--checkpoints_manifest", required=True,
                   help="Path to checkpoints_manifest.json from 02b_train_mixvae.py.")
    p.add_argument("--evaluation_results_json", required=True,
                   help="Path to evaluation_results.json from 03a_evaluate.py.")
    p.add_argument("--output_dir", required=True,
                   help="Root output directory. State files written to <output_dir>/state/.")
    p.add_argument("--kegg_toml", default="",
                   help=(
                       "Optional path to KEGG.toml. When supplied, gene-to-pathway "
                       "mapping is computed and stored in the traversal pickle. "
                       "If omitted, g_subset and pathways are stored as empty lists."
                   ))
    p.add_argument("--htree_file", default="",
                   help=(
                       "Optional path to the hierarchical taxonomy tree CSV "
                       "(e.g. tree_Mouse_ALM-VISp_2018.csv). Used to order "
                       "categories by dominant t-type taxonomy position. "
                       "Falls back to alphabetical ordering if omitted."
                   ))
    p.add_argument("--n_traversal_steps", type=int, default=50,
                   help="Number of points along the state traversal path. Default: 50.")
    p.add_argument("--traversal_std_range", type=float, default=3.0,
                   help=(
                       "Traversal spans ± (this value) × std along PC1 of the "
                       "per-category state distribution. Default: 3.0."
                   ))
    p.add_argument("--batch_size", type=int, default=5000,
                   help="Batch size for model inference. Default: 5000.")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed. Default: 0.")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# KEGG helpers
# ---------------------------------------------------------------------------

def _load_kegg(kegg_toml_path: str, gene_names: np.ndarray):
    """
    Parse KEGG.toml and return (g_subset, pathways).

    g_subset[i]  = np.array of gene indices (into gene_names) for pathway i
    pathways[i]  = pathway name string
    Only pathways with at least one gene present in gene_names are included.
    """
    import toml
    kegg = toml.load(kegg_toml_path)
    gene_set = {g: i for i, g in enumerate(gene_names)}

    g_subset, pathways = [], []
    for pathway_id, info in kegg.items():
        name = info.get("pathway_name", pathway_id).strip()
        facs_genes = info.get("facs_genes", [])
        indices = np.array([gene_set[g] for g in facs_genes if g in gene_set], dtype=int)
        if len(indices) > 0:
            g_subset.append(indices)
            pathways.append(name)

    print(f"  KEGG: {len(pathways)} pathways with genes in the model's gene list.")
    return g_subset, pathways


# ---------------------------------------------------------------------------
# Taxonomy ordering helpers
# ---------------------------------------------------------------------------

def _taxonomy_order(cluster_labels: np.ndarray,
                    predicted_label_arm0: np.ndarray,
                    nprune_indx: np.ndarray,
                    htree_file: str) -> np.ndarray:
    """
    Return category indices (0-based) ordered by their dominant t-type's
    position in the taxonomy tree.

    If htree_file is not provided or fails to load, falls back to ordering
    by alphabetically-sorted dominant t-type label.
    """
    from statistics import mode as stat_mode

    unique_ttypes = np.unique(cluster_labels)

    # For each retained category, find its dominant t-type (majority vote arm 0)
    dominant_ttype = {}
    for cc in nprune_indx:
        cells = np.where(predicted_label_arm0 == (cc + 1))[0]
        if len(cells) == 0:
            dominant_ttype[cc] = unique_ttypes[0]
        else:
            dominant_ttype[cc] = stat_mode(cluster_labels[cells])

    # Determine the ordered list of t-types
    if htree_file:
        try:
            from mmidas.utils.tree_based_analysis import get_merged_types
            _, treeobj, _ = get_merged_types(
                htree_file=htree_file,
                cells_labels=cluster_labels,
                num_classes=1,
            )
            ttype_order = [
                s.strip() for s in treeobj.child
                if (unique_ttypes == s.strip()).any()
            ]
        except Exception as e:
            print(f"  WARNING: htree_file load failed ({e}). Using alphabetical order.")
            ttype_order = sorted(unique_ttypes)
    else:
        ttype_order = sorted(unique_ttypes)

    # Rank t-types by position in ttype_order
    ttype_rank = {t: i for i, t in enumerate(ttype_order)}

    # Sort retained categories by their dominant t-type rank, then by category index
    sorted_cats = sorted(
        nprune_indx,
        key=lambda cc: (ttype_rank.get(dominant_ttype[cc], len(ttype_order)), cc),
    )
    return np.array(sorted_cats, dtype=int)


# ---------------------------------------------------------------------------
# Traversal helpers
# ---------------------------------------------------------------------------

def _traversal_for_category(
    cc: int,               # 1-indexed category id
    state_mu: np.ndarray,  # (n_arm, n_cells, state_dim)
    predicted_label: np.ndarray,  # (n_arm, n_cells) — 1-indexed
    model,                 # RNA_RNA_mixVAE (the nn.Module inside cpl_mixVAE)
    n_categories: int,
    n_arm: int,
    state_dim: int,
    n_gene: int,
    n_steps: int,
    std_range: float,
    device: torch.device,
):
    """
    Compute per-gene expression variance along a PCA-based state traversal
    for a single category cc.

    Returns
    -------
    V_g      : np.ndarray (n_arm, 1, n_gene)  — per-gene std along traversal
    V_g_std  : np.ndarray (n_arm, 1, n_gene)  — std of std (bootstrap not done; same as V_g)
    s_trav   : np.ndarray (n_arm, n_steps, state_dim) — traversal path in state space
    s_mean_arm : np.ndarray (n_arm, state_dim)  — per-arm state mean for this category
    """
    V_g     = np.zeros((n_arm, 1, n_gene))
    V_g_std = np.zeros((n_arm, 1, n_gene))
    s_trav  = np.zeros((n_arm, n_steps, state_dim))
    s_mean_arm = np.zeros((n_arm, state_dim))

    # Build category one-hot tensor (same for all arms, no grad needed)
    c_onehot = torch.zeros(1, n_categories, device=device)
    c_onehot[0, cc - 1] = 1.0  # cc is 1-indexed

    model.eval()
    with torch.no_grad():
        for arm in range(n_arm):
            cells = np.where(predicted_label[arm] == cc)[0]

            if len(cells) == 0:
                # Empty category — leave zeros
                continue

            cat_states = state_mu[arm, cells, :]  # (n_cells_c, state_dim)
            s_mean = cat_states.mean(axis=0)       # (state_dim,)
            s_mean_arm[arm] = s_mean

            # Compute traversal direction: PC1 of per-category state distribution
            if len(cells) > 1 and state_dim > 1:
                pca = PCA(n_components=1)
                pca.fit(cat_states)
                pc1_dir = pca.components_[0]   # (state_dim,)
                s_on_pc1 = pca.transform(cat_states)[:, 0]
                state_std = float(np.std(s_on_pc1)) + 1e-6
            else:
                pc1_dir = np.ones(state_dim) / np.sqrt(state_dim)
                state_std = float(np.std(cat_states) + 1e-6)

            # Traversal points in state space: mean ± std_range * std along PC1
            t_vals = np.linspace(-std_range * state_std, std_range * state_std, n_steps)
            trav_points = s_mean[None, :] + t_vals[:, None] * pc1_dir[None, :]
            # shape: (n_steps, state_dim)
            s_trav[arm] = trav_points

            # Decode each traversal point
            recon_trav = np.zeros((n_steps, n_gene))
            for t_idx in range(n_steps):
                s_t = torch.FloatTensor(trav_points[t_idx : t_idx + 1]).to(device)
                recon = model.decoder(c_onehot, s_t, arm)
                recon_trav[t_idx] = recon.cpu().numpy().flatten()

            V_g[arm, 0, :]     = np.std(recon_trav, axis=0)
            V_g_std[arm, 0, :] = V_g[arm, 0, :]   # placeholder — no bootstrap

    return V_g, V_g_std, s_trav, s_mean_arm


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)
    np.random.seed(args.seed)

    state_dir = os.path.join(args.output_dir, "state")
    os.makedirs(state_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Load manifests
    # ------------------------------------------------------------------
    with open(args.checkpoints_manifest) as fh:
        ckpt_manifest = json.load(fh)
    with open(args.evaluation_results_json) as fh:
        eval_results = json.load(fh)

    n_categories  = ckpt_manifest["n_categories"]
    state_dim     = ckpt_manifest["state_dim"]
    n_arm         = ckpt_manifest["n_arm"]
    latent_dim    = ckpt_manifest["latent_dim"]
    fc_dim        = ckpt_manifest.get("fc_dim", 100)
    x_drop        = ckpt_manifest.get("x_drop", 0.0)
    s_drop        = ckpt_manifest.get("s_drop", 0.0)
    temp          = ckpt_manifest.get("temp", 1.0)
    tau           = ckpt_manifest.get("tau", 0.1)
    beta          = ckpt_manifest.get("beta", 1.0)
    lam           = ckpt_manifest.get("lam", 1.0)
    lam_pc        = ckpt_manifest.get("lam_pc", 1.0)
    hard          = ckpt_manifest.get("hard", False)
    variational   = ckpt_manifest.get("variational", True)
    training_mode = ckpt_manifest.get("training_mode", "MSE")
    n_gene        = ckpt_manifest["n_gene"]

    model_order    = eval_results["model_order"]
    selected_model = eval_results["selected_model"]

    print(f"model_order={model_order}, n_arm={n_arm}, state_dim={state_dim}, n_gene={n_gene}")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"\nLoading data: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x            = data["log1p"]
    cluster_labels = data["cluster"]
    gene_names   = data["gene_id"]   # array of gene symbols
    n_cells      = x.shape[0]
    print(f"  {n_cells} cells, {len(gene_names)} genes")

    # ------------------------------------------------------------------
    # KEGG gene-to-pathway mapping (optional)
    # ------------------------------------------------------------------
    if args.kegg_toml and os.path.exists(args.kegg_toml):
        print(f"\nLoading KEGG from: {args.kegg_toml}")
        g_subset, pathways = _load_kegg(args.kegg_toml, gene_names)
    else:
        print("\nNo KEGG file provided — g_subset and pathways will be empty.")
        g_subset, pathways = [], []

    # ------------------------------------------------------------------
    # Model inference
    # ------------------------------------------------------------------
    print("\nRunning model inference ...")
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
    mixvae.variational = False

    all_loader = DataLoader(
        TensorDataset(
            torch.FloatTensor(x),
            torch.FloatTensor(range(n_cells)),
        ),
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )

    outcome = summarize_inference(
        cpl_mixVAE=mixvae,
        files=[selected_model],
        data=all_loader,
        saving_folder="",
    )

    # outcome keys used here:
    # state_mu:      (n_arm, n_cells, state_dim)
    # pred_label:    list[checkpoint] → [0] = (n_arm, n_cells) 1-indexed
    # nprune_indx:   array of retained 0-based category indices
    state_mu        = outcome["state_mu"]           # (n_arm, n_cells, state_dim)
    predicted_label = outcome["pred_label"][0]      # (n_arm, n_cells) 1-indexed
    nprune_indx     = outcome["nprune_indx"]        # 0-based retained indices

    # Active categories: 1-indexed
    c_cat = np.array([cc + 1 for cc in nprune_indx], dtype=int)
    n_cats = len(c_cat)
    print(f"  {n_cats} active categories to process.")

    # ------------------------------------------------------------------
    # Taxonomy ordering
    # ------------------------------------------------------------------
    print("\nComputing taxonomy ordering ...")
    taxonomy_order = _taxonomy_order(
        cluster_labels=np.array([c.strip() for c in cluster_labels]),
        predicted_label_arm0=predicted_label[0],
        nprune_indx=nprune_indx,
        htree_file=args.htree_file,
    )
    taxonomy_path = os.path.join(
        args.output_dir, f"taxonomy_order_K_{model_order}.npy"
    )
    np.save(taxonomy_path, taxonomy_order)
    print(f"  Saved: {taxonomy_path}")

    # ------------------------------------------------------------------
    # State traversal computation
    # ------------------------------------------------------------------
    print(f"\nComputing state traversals ({args.n_traversal_steps} steps, "
          f"±{args.traversal_std_range}σ along PC1) ...")

    device = mixvae.device
    model  = mixvae.model   # the nn.Module (RNA_RNA_mixVAE)

    V_g_mean = np.zeros((n_cats, n_arm, 1, n_gene))
    V_g_std  = np.zeros((n_cats, n_arm, 1, n_gene))
    s_travers = []      # list of (n_arm, n_steps, state_dim) per category
    s_mu_list = []      # list of (n_arm, state_dim) per category

    for ic, cc in enumerate(c_cat):
        V_g, Vg_s, s_t, s_m = _traversal_for_category(
            cc=int(cc),
            state_mu=state_mu,
            predicted_label=predicted_label,
            model=model,
            n_categories=n_categories,
            n_arm=n_arm,
            state_dim=state_dim,
            n_gene=n_gene,
            n_steps=args.n_traversal_steps,
            std_range=args.traversal_std_range,
            device=device,
        )
        V_g_mean[ic] = V_g
        V_g_std[ic]  = Vg_s
        s_travers.append(s_t)
        s_mu_list.append(s_m)

        if (ic + 1) % max(1, n_cats // 10) == 0 or ic == n_cats - 1:
            print(f"  {ic + 1}/{n_cats} categories done.")

    # ------------------------------------------------------------------
    # Save traversal_pca pickle
    # ------------------------------------------------------------------
    traversal_path = os.path.join(
        state_dir, f"traversal_pca_K_{model_order}.pickle"
    )
    with open(traversal_path, "wb") as fh:
        pickle.dump(
            {
                "V_g_mean": V_g_mean,
                "V_g_std":  V_g_std,
                "g_subset": g_subset,
                "c_cat":    c_cat,
                "pathways": pathways,
            },
            fh,
            protocol=4,
        )
    print(f"\n  Saved: {traversal_path}")

    # ------------------------------------------------------------------
    # Save state_mu_pca pickle
    # ------------------------------------------------------------------
    state_mu_path = os.path.join(
        state_dir, f"state_mu_pca_K_{model_order}.pickle"
    )
    with open(state_mu_path, "wb") as fh:
        pickle.dump(
            {
                "s_mu":     np.array(s_mu_list),  # (n_cats, n_arm, state_dim)
                "s_travers": s_travers,            # list of (n_arm, n_steps, state_dim)
                "c_cat":    c_cat,
            },
            fh,
            protocol=4,
        )
    print(f"  Saved: {state_mu_path}")

    # ------------------------------------------------------------------
    # Manifest
    # ------------------------------------------------------------------
    manifest = {
        "output_dir":       os.path.abspath(args.output_dir),
        "state_dir":        os.path.abspath(state_dir),
        "model_order":      model_order,
        "n_categories":     n_categories,
        "n_arm":            n_arm,
        "state_dim":        state_dim,
        "n_gene":           n_gene,
        "n_cats":           n_cats,
        "n_traversal_steps": args.n_traversal_steps,
        "traversal_pickle": traversal_path,
        "state_mu_pickle":  state_mu_path,
        "taxonomy_order":   taxonomy_path,
        "n_pathways":       len(pathways),
    }
    manifest_path = os.path.join(args.output_dir, "traversal_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nTraversal prep complete.")
    print(f"  {n_cats} categories, {len(pathways)} KEGG pathways")
    print(f"  Manifest: {manifest_path}")


if __name__ == "__main__":
    sys.exit(main())
