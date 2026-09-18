#!/usr/bin/env python3
"""
03b_classify.py — MMIDAS Clusterability: Classification & Silhouette Analysis

Measures the post-hoc identifiability of cell-type labels through K-fold
Random Forest classification and silhouette scoring. Produces pickle files
consumed by 04_clusterability.py.

Two embedding types are compared:
  - Linear:     PCA of raw log-CPM expression (n_pca components)
  - Non-linear: Per-arm MMIDAS low-dimensional embedding (latent_dim dims)

Three label sets are classified:
  - Ttype:      Original cluster labels from obs['cluster']
  - ConsType:   MMIDAS-inferred category labels (per arm, from inference)

Output pickle files (written to <output_dir>/clustering/):
  Global (no arm suffix — PCA features, averaged over arms for MMIDAS):
    Ttype_classification_K_{n_ttype}_nFeature_{n_pca}_{date}.p
  Per arm:
    Ttype_classification_K_{n_ttype}_nFeature_{latent_dim}_arm_{arm}_{date}.p
    ConsType_classification_K_{model_order}_nFeature_{n_pca}_arm_{arm}_{date}.p
    ConsType_classification_K_{model_order}_nFeature_{latent_dim}_arm_{arm}_{date}.p

Each pickle contains:
  acc_T_adj  — np.array of shape (k_fold,): per-fold RF accuracy
  sc_T       — list of np.arrays: per-cluster mean silhouette scores
  conf_mat   — np.array: aggregated confusion matrix across folds

Usage:
  python3 03b_classify.py \\
    --anndata_h5ad            Mouse_ALM-VISp_cpm.h5ad \\
    --checkpoints_manifest    results/run_1/checkpoints_manifest.json \\
    --evaluation_results_json results/run_1/evaluation_results.json \\
    --output_dir              results/run_1

For a fast smoke test use: --k_fold 2 --n_pca 5
"""

import argparse
import json
import os
import pickle
import sys
from datetime import datetime

import matplotlib
matplotlib.use("Agg")

import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.metrics import silhouette_samples
from sklearn.model_selection import KFold

from torch.utils.data import DataLoader, TensorDataset
import torch

from mmidas.cpl_mixvae import cpl_mixVAE
from mmidas.eval import summarize_inference
from mmidas.utils.data_tools import load_data


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Run RF classification and silhouette analysis on MMIDAS embeddings."
    )

    p.add_argument("--anndata_h5ad", required=True,
                   help="Path to the .h5ad file (output of 01_data_prep.py).")
    p.add_argument("--checkpoints_manifest", required=True,
                   help="Path to checkpoints_manifest.json from 02b_train_mixvae.py.")
    p.add_argument("--evaluation_results_json", required=True,
                   help="Path to evaluation_results.json from 03a_evaluate.py.")
    p.add_argument("--output_dir", required=True,
                   help="Root output directory. Pickles are written to <output_dir>/clustering/.")
    p.add_argument("--n_pca", type=int, default=100,
                   help="Number of PCA components for the linear embedding. Default: 100.")
    p.add_argument("--k_fold", type=int, default=10,
                   help="Number of K-fold CV splits for RF classification. Default: 10.")
    p.add_argument("--batch_size", type=int, default=5000,
                   help="Batch size for model inference. Default: 5000.")
    p.add_argument("--seed", type=int, default=0,
                   help="Random seed for KFold and RF. Default: 0.")
    p.add_argument("--date_tag", default="",
                   help=(
                       "Date string embedded in output filenames (e.g. 20231121). "
                       "Defaults to today's date (YYYYMMDD)."
                   ))

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Classification helpers
# ---------------------------------------------------------------------------

def _rf_classify(features: np.ndarray, labels: np.ndarray,
                 k_fold: int, seed: int):
    """
    K-fold Random Forest classification.

    Returns
    -------
    acc      : np.ndarray, shape (k_fold,)
    conf_mat : np.ndarray, shape (n_classes, n_classes) — summed across folds
    """
    kf = KFold(n_splits=k_fold, shuffle=True, random_state=seed)
    classes = np.unique(labels)
    n_cls = len(classes)
    acc = []
    conf_mat = np.zeros((n_cls, n_cls), dtype=int)

    for train_idx, test_idx in kf.split(features):
        clf = RandomForestClassifier(random_state=seed)
        clf.fit(features[train_idx], labels[train_idx])
        pred = clf.predict(features[test_idx])
        acc.append(accuracy_score(labels[test_idx], pred))
        conf_mat += confusion_matrix(labels[test_idx], pred, labels=classes)

    return np.array(acc), conf_mat


def _silhouette_per_cluster(features: np.ndarray, labels: np.ndarray):
    """
    Per-cluster mean silhouette scores.

    Returns list with one element (to match the notebook's sc_T[-1] convention).
    Returns [np.zeros(n_classes)] if there is only one unique label
    (silhouette is undefined).
    """
    classes = np.unique(labels)
    if len(classes) < 2:
        return [np.zeros(len(classes))]

    sample_scores = silhouette_samples(features, labels)
    mean_per_cluster = np.array([
        np.mean(sample_scores[labels == c]) for c in classes
    ])
    return [mean_per_cluster]


def _save_pickle(path: str, acc: np.ndarray, sc_T: list, conf_mat: np.ndarray):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as fh:
        pickle.dump({"acc_T_adj": acc, "sc_T": sc_T, "conf_mat": conf_mat}, fh, protocol=4)
    print(f"  Saved: {os.path.basename(path)}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    date_tag = args.date_tag or datetime.now().strftime("%Y%m%d")
    clustering_dir = os.path.join(args.output_dir, "clustering")
    os.makedirs(clustering_dir, exist_ok=True)

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

    print(f"model_order={model_order}, n_arm={n_arm}, latent_dim={latent_dim}")

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"\nLoading data: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x = data["log1p"]                           # (n_cells, n_gene)
    cluster_labels = data["cluster"]            # obs['cluster'] string labels
    n_cells = x.shape[0]
    n_ttype = len(np.unique(cluster_labels))
    print(f"  {n_cells} cells, {n_ttype} unique cluster labels")

    # ------------------------------------------------------------------
    # Linear embedding: PCA
    # ------------------------------------------------------------------
    n_pca = min(args.n_pca, n_cells - 1, n_gene)
    print(f"\nComputing PCA ({n_pca} components)...")
    pca = PCA(n_components=n_pca, random_state=args.seed)
    x_pca = pca.fit_transform(x)
    print(f"  PCA shape: {x_pca.shape}")

    # ------------------------------------------------------------------
    # Non-linear embedding: MMIDAS per-arm lowD representation
    # ------------------------------------------------------------------
    print(f"\nRunning model inference to get per-arm embeddings...")
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

    all_tensor = torch.FloatTensor(x)
    all_idx    = torch.FloatTensor(range(n_cells))
    all_loader = DataLoader(
        TensorDataset(all_tensor, all_idx),
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

    # outcome['lowD_x']    shape: (n_arm, n_cells, latent_dim)
    # outcome['pred_label'] is a list indexed by checkpoint; [0] = selected model
    #                        each entry shape: (n_arm, n_cells) — 1-indexed categories
    lowD_x          = outcome["lowD_x"]            # (n_arm, n_cells, latent_dim)
    predicted_label = outcome["pred_label"][0]     # (n_arm, n_cells)

    # ------------------------------------------------------------------
    # Label sets
    # ------------------------------------------------------------------
    # Integer-encode cluster labels for RF (RF works on string labels too,
    # but integer encoding makes confusion_matrix consistent).
    unique_ttypes = np.unique(cluster_labels)
    ttype_int = np.array([
        np.where(unique_ttypes == c)[0][0] for c in cluster_labels
    ])

    produced_files = []

    # ------------------------------------------------------------------
    # 1. Ttype × PCA  (global, no arm suffix)
    # ------------------------------------------------------------------
    print(f"\n[1/4] Ttype × PCA (K={args.k_fold} folds) ...")
    acc, conf = _rf_classify(x_pca, ttype_int, args.k_fold, args.seed)
    sc = _silhouette_per_cluster(x_pca, ttype_int)
    fname = os.path.join(
        clustering_dir,
        f"Ttype_classification_K_{n_ttype}_nFeature_{n_pca}_{date_tag}.p"
    )
    _save_pickle(fname, acc, sc, conf)
    produced_files.append(fname)

    for arm in range(n_arm):
        arm_emb = lowD_x[arm]   # (n_cells, latent_dim)
        cons_int = (predicted_label[arm] - 1).astype(int)  # 0-indexed

        # ------------------------------------------------------------------
        # 2. Ttype × per-arm MMIDAS embedding
        # ------------------------------------------------------------------
        print(f"\n[2/4] Ttype × MMIDAS arm {arm} ...")
        acc, conf = _rf_classify(arm_emb, ttype_int, args.k_fold, args.seed)
        sc = _silhouette_per_cluster(arm_emb, ttype_int)
        fname = os.path.join(
            clustering_dir,
            f"Ttype_classification_K_{n_ttype}_nFeature_{latent_dim}_arm_{arm}_{date_tag}.p"
        )
        _save_pickle(fname, acc, sc, conf)
        produced_files.append(fname)

        # ------------------------------------------------------------------
        # 3. ConsType × PCA  (per arm)
        # ------------------------------------------------------------------
        print(f"\n[3/4] ConsType × PCA arm {arm} ...")
        acc, conf = _rf_classify(x_pca, cons_int, args.k_fold, args.seed)
        sc = _silhouette_per_cluster(x_pca, cons_int)
        fname = os.path.join(
            clustering_dir,
            f"ConsType_classification_K_{model_order}_nFeature_{n_pca}_arm_{arm}_{date_tag}.p"
        )
        _save_pickle(fname, acc, sc, conf)
        produced_files.append(fname)

        # ------------------------------------------------------------------
        # 4. ConsType × per-arm MMIDAS embedding
        # ------------------------------------------------------------------
        print(f"\n[4/4] ConsType × MMIDAS arm {arm} ...")
        acc, conf = _rf_classify(arm_emb, cons_int, args.k_fold, args.seed)
        sc = _silhouette_per_cluster(arm_emb, cons_int)
        fname = os.path.join(
            clustering_dir,
            f"ConsType_classification_K_{model_order}_nFeature_{latent_dim}_arm_{arm}_{date_tag}.p"
        )
        _save_pickle(fname, acc, sc, conf)
        produced_files.append(fname)

    # ------------------------------------------------------------------
    # Manifest
    # ------------------------------------------------------------------
    manifest = {
        "output_dir":    os.path.abspath(args.output_dir),
        "clustering_dir": os.path.abspath(clustering_dir),
        "model_order":   model_order,
        "n_ttype":       n_ttype,
        "n_pca":         n_pca,
        "latent_dim":    latent_dim,
        "n_arm":         n_arm,
        "k_fold":        args.k_fold,
        "date_tag":      date_tag,
        "pickles":       produced_files,
    }
    manifest_path = os.path.join(args.output_dir, "classify_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nClassification complete.")
    print(f"  {len(produced_files)} pickle(s) written to {clustering_dir}")
    print(f"  Manifest: {manifest_path}")


if __name__ == "__main__":
    sys.exit(main())
