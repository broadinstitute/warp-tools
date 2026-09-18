#!/usr/bin/env python3
"""
02a_train_augmenter.py — MMIDAS UDAGAN Augmenter Training

Trains the VAE-GAN data augmenter (UDAGAN) on a prepared AnnData .h5ad file.
The resulting checkpoint can be passed to 02b_train_mixvae.py via
--augmenter_checkpoint to improve cpl-mixVAE training on small datasets.

Steps:
  1. Load the .h5ad file via mmidas.utils.data_tools.load_data().
  2. Build a binary-masked DataLoader via mmidas.utils.augmentation.get_loader().
  3. Initialise vae_gan and run training for n_epoch epochs.
  4. Write a JSON manifest with the checkpoint path and hyperparameters.

Usage:
  python3 02a_train_augmenter.py \\
    --anndata_h5ad  Mouse_ALM-VISp_cpm.h5ad \\
    --output_dir    results/augmenter \\
    [--cuda]

For a fast smoke test use: --n_epoch 2 --batch_size 16
"""

import argparse
import json
import os
import sys
from pathlib import Path

# Set headless matplotlib backend BEFORE any mmidas/torch imports.
# vae_gan.train() calls plt.savefig() internally.
import matplotlib
matplotlib.use("Agg")

import torch
from mmidas.vaegan import vae_gan
from mmidas.utils.data_tools import load_data
from mmidas.utils.augmentation import get_loader


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Train the MMIDAS UDAGAN VAE-GAN data augmenter."
    )

    # --- I/O ------------------------------------------------------------------
    p.add_argument(
        "--anndata_h5ad", required=True,
        help="Path to the input .h5ad file (output of 01_data_prep.py)."
    )
    p.add_argument(
        "--output_dir", required=True,
        help=(
            "Directory for all outputs. Checkpoint and loss_curve.png are "
            "written here. Created if it does not exist."
        )
    )
    p.add_argument(
        "--tag", default="mmidas",
        help=(
            "Short string embedded in the checkpoint filename: "
            "RNA_augmenter_<tag>_<timestamp>.pth. Default: 'mmidas'."
        )
    )

    # --- Architecture ---------------------------------------------------------
    p.add_argument(
        "--z_dim", type=int, default=10,
        help="Latent space dimension of the augmenter encoder. Default: 10."
    )
    p.add_argument(
        "--noise_dim", type=int, default=50,
        help="Additive noise dimension. Default: 50."
    )
    p.add_argument(
        "--fc_dim", type=int, default=500,
        help="Hidden layer width. Default: 500."
    )
    p.add_argument(
        "--x_drop", type=float, default=0.25,
        help="Input dropout probability. Default: 0.25."
    )
    p.add_argument(
        "--affine", action="store_true",
        help="Use affine transformation in batch normalisation layers."
    )
    p.add_argument(
        "--momentum", type=float, default=0.01,
        help="Batch normalisation momentum. Default: 0.01."
    )

    # --- Training -------------------------------------------------------------
    p.add_argument(
        "--n_epoch", type=int, default=10,
        help="Number of training epochs. Default: 10."
    )
    p.add_argument(
        "--lr", type=float, default=1e-3,
        help="Adam optimizer learning rate. Default: 0.001."
    )
    p.add_argument(
        "--alpha", type=float, default=0.2,
        help="Triplet loss margin. Default: 0.2."
    )
    p.add_argument(
        "--lam", type=float, nargs=4, default=[1.0, 0.5, 0.1, 0.5],
        metavar=("LAM_GEN", "LAM_TRIP", "LAM_Z", "LAM_RECON"),
        help=(
            "Four loss weights: [generator, triplet, z-distance, reconstruction]. "
            "Default: 1.0 0.5 0.1 0.5."
        )
    )
    p.add_argument(
        "--batch_size", type=int, default=512,
        help="Mini-batch size. Default: 512."
    )
    p.add_argument(
        "--cuda", action="store_true",
        help="Use GPU (CUDA) if available. Falls back to CPU with a warning."
    )

    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None):
    args = parse_args(argv)

    # ------------------------------------------------------------------
    # Device selection
    # ------------------------------------------------------------------
    if args.cuda:
        if not torch.cuda.is_available():
            print(
                "WARNING: --cuda requested but CUDA is not available. "
                "Falling back to CPU."
            )
            device = torch.device("cpu")
        else:
            free_gpus = [
                i for i in range(torch.cuda.device_count())
                if (
                    torch.cuda.get_device_properties(i).total_memory
                    - torch.cuda.memory_allocated(i) > 0
                )
            ]
            if not free_gpus:
                raise RuntimeError("--cuda requested but no free CUDA devices found.")
            device = torch.device(f"cuda:{free_gpus[0]}")
            print(f"Using GPU: {torch.cuda.get_device_name(device)}")
    else:
        device = torch.device("cpu")
        print("Using CPU.")

    # ------------------------------------------------------------------
    # Output directory
    # vae_gan.train() uses Path division (self.folder / filename) so
    # saving_folder must be a pathlib.Path, not a plain string.
    # ------------------------------------------------------------------
    saving_folder = Path(args.output_dir)
    saving_folder.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"Loading data from: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x = data["log1p"]   # dense numpy array: n_cells × n_genes
    n_gene = x.shape[1]
    print(f"  {x.shape[0]} cells × {n_gene} genes loaded.")

    # ------------------------------------------------------------------
    # DataLoader
    # ------------------------------------------------------------------
    print("Building augmenter DataLoader...")
    data_loader = get_loader(x=x, batch_size=args.batch_size, training=True)

    # ------------------------------------------------------------------
    # Model initialisation
    # ------------------------------------------------------------------
    augmenter = vae_gan(saving_folder=saving_folder, device=device)
    augmenter.init_model(
        input_dim=n_gene,
        z_dim=args.z_dim,
        noise_dim=args.noise_dim,
        fc_dim=args.fc_dim,
        x_drop=args.x_drop,
        affine=args.affine,
        momentum=args.momentum,
    )

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------
    print(
        f"\nStarting augmenter training: "
        f"n_epoch={args.n_epoch}, batch_size={args.batch_size}, "
        f"z_dim={args.z_dim}, noise_dim={args.noise_dim}"
    )

    checkpoint_filename = augmenter.train(
        dataloader=data_loader,
        n_epoch=args.n_epoch,
        lr=args.lr,
        alpha=args.alpha,
        lam=args.lam,
        tag=args.tag,
    )

    checkpoint_path = str(saving_folder / checkpoint_filename)

    # ------------------------------------------------------------------
    # Manifest
    # ------------------------------------------------------------------
    manifest = {
        "output_dir":         str(saving_folder.resolve()),
        "checkpoint":         checkpoint_path,
        "n_gene":             n_gene,
        "z_dim":              args.z_dim,
        "noise_dim":          args.noise_dim,
        "fc_dim":             args.fc_dim,
    }
    manifest_path = saving_folder / "augmenter_manifest.json"
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nAugmenter training complete.")
    print(f"  Checkpoint: {checkpoint_path}")
    print(f"  Manifest:   {manifest_path}")


if __name__ == "__main__":
    sys.exit(main())
