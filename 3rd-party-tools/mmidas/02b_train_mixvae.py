#!/usr/bin/env python3
"""
02b_train_mixvae.py — MMIDAS cpl-mixVAE Training

Trains the Coupled Mixture VAE (cpl-mixVAE) model on a prepared AnnData .h5ad
file produced by 01_data_prep.py.

Steps:
  1. Load the .h5ad file via mmidas.utils.data_tools.load_data().
  2. Build train/test DataLoaders.
  3. Instantiate cpl_mixVAE and run initial training for n_epoch epochs.
  4. Run iterative pruning for up to max_prun_it rounds (n_epoch_p epochs each),
     eliminating the lowest-consensus category each round.
  5. Write a JSON manifest listing all produced .pth checkpoint files.

Usage:
  python3 02b_train_mixvae.py \\
    --anndata_h5ad  Mouse_ALM-VISp_cpm.h5ad \\
    --output_dir    results/run_1 \\
    [--augmenter_checkpoint augmenter.pth] \\
    [--cuda]

All other arguments are optional hyperparameters with production defaults.
For a fast smoke test use: --n_epoch 2 --n_epoch_p 1 --max_prun_it 1 --batch_size 16
"""

import argparse
import json
import os
import sys

# Set headless matplotlib backend BEFORE any mmidas/torch imports.
# cpl_mixvae.py calls plt.savefig() internally during training, so this must
# be set before the module is imported.
import matplotlib
matplotlib.use("Agg")

import torch
from mmidas.cpl_mixvae import cpl_mixVAE
from mmidas.utils.data_tools import load_data, get_loaders


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Train the MMIDAS cpl-mixVAE model on a prepared .h5ad file."
    )

    # --- I/O ------------------------------------------------------------------
    p.add_argument(
        "--anndata_h5ad", required=True,
        help="Path to the input .h5ad file (output of 01_data_prep.py)."
    )
    p.add_argument(
        "--output_dir", required=True,
        help=(
            "Root directory for all outputs. Checkpoints are written to "
            "<output_dir>/model/; figures go to <output_dir>/. "
            "Created if it does not exist."
        )
    )
    p.add_argument(
        "--augmenter_checkpoint", default="",
        help=(
            "Optional path to a pre-trained UDAGAN augmenter .pth checkpoint "
            "(produced by 02a_train_augmenter.py). If omitted, no augmentation "
            "is applied."
        )
    )

    # --- Architecture ---------------------------------------------------------
    p.add_argument(
        "--n_categories", type=int, default=120,
        help="Initial (upper-bound) number of cell-type categories. Default: 120."
    )
    p.add_argument(
        "--state_dim", type=int, default=2,
        help="Dimension of the continuous state latent variable. Default: 2."
    )
    p.add_argument(
        "--n_arm", type=int, default=2,
        help="Number of coupled encoder arms. Default: 2."
    )
    p.add_argument(
        "--latent_dim", type=int, default=10,
        help="Dimension of the low-dimensional embedding (lowD_dim). Default: 10."
    )
    p.add_argument(
        "--fc_dim", type=int, default=100,
        help="Hidden layer width. Default: 100."
    )
    p.add_argument(
        "--x_drop", type=float, default=0.25,
        help="Input-layer dropout probability. Default: 0.25."
    )
    p.add_argument(
        "--s_drop", type=float, default=0.0,
        help="State-variable dropout probability. Default: 0.0."
    )
    p.add_argument(
        "--temp", type=float, default=1.0,
        help="Gumbel-softmax temperature. Default: 1.0."
    )
    p.add_argument(
        "--tau", type=float, default=0.1,
        help="Softmax temperature for the categorical variable. Default: 0.1."
    )
    p.add_argument(
        "--beta", type=float, default=1.0,
        help="KL divergence regularization coefficient. Default: 1.0."
    )
    p.add_argument(
        "--lam", type=float, default=1.0,
        help="Coupling factor between arms. Default: 1.0."
    )
    p.add_argument(
        "--lam_pc", type=float, default=1.0,
        help="Coupling factor for the reference arm. Default: 1.0."
    )
    p.add_argument(
        "--hard", action="store_true",
        help="Use hard (one-hot) Gumbel-softmax encoding instead of soft."
    )
    p.add_argument(
        "--variational", type=lambda x: x.lower() != "false", default=True,
        help=(
            "Enable variational mode (True/False). Default: True. "
            "Pass --variational False to disable."
        )
    )
    p.add_argument(
        "--training_mode", default="MSE", choices=["MSE", "ZINB"],
        help="Reconstruction loss function. Default: MSE."
    )

    # --- Training -------------------------------------------------------------
    p.add_argument(
        "--n_epoch", type=int, default=10000,
        help="Training epochs before pruning begins. Default: 10000."
    )
    p.add_argument(
        "--n_epoch_p", type=int, default=1000,
        help="Training epochs per pruning iteration. Default: 1000."
    )
    p.add_argument(
        "--min_con", type=float, default=0.99,
        help=(
            "Minimum inter-arm consensus to retain a category. "
            "Pruning stops early if all remaining categories exceed this. Default: 0.99."
        )
    )
    p.add_argument(
        "--max_prun_it", type=int, default=14,
        help="Maximum number of pruning iterations. Default: 14."
    )
    p.add_argument(
        "--lr", type=float, default=1e-3,
        help="Adam optimizer learning rate. Default: 0.001."
    )
    p.add_argument(
        "--batch_size", type=int, default=5000,
        help="Mini-batch size. Default: 5000."
    )
    p.add_argument(
        "--n_aug_smp", type=int, default=0,
        help=(
            "Number of augmented copies appended per training sample. "
            "Only used when --augmenter_checkpoint is provided. Default: 0."
        )
    )
    p.add_argument(
        "--train_size", type=float, default=0.9,
        help="Fraction of data used for training; remainder is test. Default: 0.9."
    )
    p.add_argument(
        "--seed", type=int, default=0,
        help="Random seed for the train/test split. Default: 0."
    )
    p.add_argument(
        "--cuda", action="store_true",
        help="Use GPU (CUDA) if available. Falls back to CPU with a warning if not."
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
    # ------------------------------------------------------------------
    model_dir = os.path.join(args.output_dir, "model")
    os.makedirs(model_dir, exist_ok=True)

    # ------------------------------------------------------------------
    # Optional augmenter
    # ------------------------------------------------------------------
    augmenter = []
    if args.augmenter_checkpoint:
        from mmidas.vaegan import vae_gan  # imported lazily — not needed for CPU runs
        print(f"Loading augmenter from: {args.augmenter_checkpoint}")
        aug_dir = os.path.dirname(os.path.abspath(args.augmenter_checkpoint))
        aug_vaegan = vae_gan(saving_folder=aug_dir, device=device)
        aug_vaegan.load_model(args.augmenter_checkpoint)
        augmenter = aug_vaegan.netA

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    print(f"Loading data from: {args.anndata_h5ad}")
    data = load_data(file=args.anndata_h5ad)
    x = data["log1p"]   # numpy array: n_cells × n_genes (dense, from .toarray())
    n_gene = x.shape[1]
    print(f"  {x.shape[0]} cells × {n_gene} genes loaded.")

    # ------------------------------------------------------------------
    # Data loaders
    # ------------------------------------------------------------------
    print("Building data loaders...")
    # get_loaders returns: alldata_loader, train_loader, test_loader,
    #                      validation_loader, test_ind, val_ind
    _, train_loader, test_loader, _, _, _ = get_loaders(
        x=x,
        batch_size=args.batch_size,
        n_aug_smp=args.n_aug_smp,
        netA=augmenter if augmenter else None,
        additional_val=False,
        seed=args.seed,
    )
    print(
        f"  Train batches: {len(train_loader)}, "
        f"Test samples: {len(test_loader.dataset)}"
    )

    # ------------------------------------------------------------------
    # Model initialisation
    # ------------------------------------------------------------------
    mixvae = cpl_mixVAE(
        saving_folder=args.output_dir,
        augmenter=augmenter,
        device=device,
    )

    mixvae.init_model(
        n_categories=args.n_categories,
        state_dim=args.state_dim,
        input_dim=n_gene,
        fc_dim=args.fc_dim,
        lowD_dim=args.latent_dim,
        x_drop=args.x_drop,
        s_drop=args.s_drop,
        n_arm=args.n_arm,
        temp=args.temp,
        hard=args.hard,
        variational=args.variational,
        tau=args.tau,
        lam=args.lam,
        lam_pc=args.lam_pc,
        beta=args.beta,
        ref_prior=False,
        mode=args.training_mode,
        trained_model="",   # fresh training run (no warm-start)
        n_pr=0,
    )

    # ------------------------------------------------------------------
    # Training
    # ------------------------------------------------------------------
    print(
        f"\nStarting training: "
        f"n_epoch={args.n_epoch}, "
        f"n_epoch_p={args.n_epoch_p}, "
        f"max_prun_it={args.max_prun_it}, "
        f"n_categories={args.n_categories}"
    )

    mixvae.train(
        train_loader=train_loader,
        test_loader=test_loader,
        n_epoch=args.n_epoch,
        n_epoch_p=args.n_epoch_p,
        lr=args.lr,
        min_con=args.min_con,
        max_prun_it=args.max_prun_it,
    )

    # ------------------------------------------------------------------
    # Checkpoint manifest
    # ------------------------------------------------------------------
    # Collect all .pth files written by cpl_mixvae.train() into model_dir.
    # This JSON file is the stable contract between this script and the
    # downstream 03a_evaluate.py task.
    checkpoints = sorted([
        os.path.join(model_dir, f)
        for f in os.listdir(model_dir)
        if f.endswith(".pth")
    ])

    manifest = {
        # I/O
        "output_dir":    os.path.abspath(args.output_dir),
        "checkpoints":   checkpoints,
        # Architecture — all fields required to reconstruct the model
        # exactly when loading checkpoints in 03a_evaluate.py
        "n_categories":  args.n_categories,
        "state_dim":     args.state_dim,
        "n_arm":         args.n_arm,
        "latent_dim":    args.latent_dim,
        "fc_dim":        args.fc_dim,
        "x_drop":        args.x_drop,
        "s_drop":        args.s_drop,
        "temp":          args.temp,
        "tau":           args.tau,
        "beta":          args.beta,
        "lam":           args.lam,
        "lam_pc":        args.lam_pc,
        "hard":          args.hard,
        "variational":   args.variational,
        "training_mode": args.training_mode,
        "n_gene":        n_gene,
    }
    manifest_path = os.path.join(args.output_dir, "checkpoints_manifest.json")
    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(f"\nTraining complete.")
    print(f"  {len(checkpoints)} checkpoint(s):")
    for ckpt in checkpoints:
        print(f"    {ckpt}")
    print(f"  Manifest: {manifest_path}")


if __name__ == "__main__":
    sys.exit(main())
