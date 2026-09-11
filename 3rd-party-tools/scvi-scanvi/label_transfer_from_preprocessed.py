"""scANVI label transfer from PREPROCESSED inputs — the post-preprocessing half of the pipeline.

PreprocessFilter (in scANVI.wdl) has already filtered/batch-labeled/modality-tagged the GEX,
reference, and (optional) ATAC-activity h5ads. This script runs everything after that: obtain the
SCANVI model (LOAD a supplied .tar.gz, else TRAIN SCVI+SCANVI), transfer labels, optionally record
per-cell max posterior, propagate labels back to the matrices, finalize metadata, and write outputs.

Kept separate from main() (which re-runs preprocessing) so the WDL task is a single command instead
of a large inline heredoc, and so the identical command runs unchanged on GPU or CPU nodes —
scvi-tools auto-detects the accelerator (accelerator="auto"), so no device handling is needed here.
Living in the container with run_multi_model/run_gex_only_model also means the training kwargs
(max_epochs, batch_size) can never skew against the WDL.

Writes to the working directory:
  <input_id>_SCANVI_predictions.h5ad
  <input_id>_gex_annotated_matrix.h5ad
  <input_id>_atac_annotated_matrix.h5ad   (multiome mode only)
  <input_id>_scanvi_model.tar.gz          (minimal SCANVI model, no bundled adata)

  python3 label_transfer_from_preprocessed.py --gex g.h5ad --ref r.h5ad --input-id ID \
      [--atac a.h5ad] [--scanvi-model m.tar.gz] [--max-epochs N] [--batch-size N] \
      [--output-max-probability]
"""
import argparse
import glob
import os
import tarfile
import time

import anndata as ad
import scanpy as sc

# PyTorch 2.6 flipped torch.load(weights_only) to True by default, which breaks scvi-tools 1.2 model
# loading (the checkpoint pickles numpy globals). Any model loaded here is produced by this same
# pipeline (trusted), so force weights_only=False for all torch.load calls.
import torch
_orig_torch_load = torch.load
def _torch_load_compat(*a, **k):
    k.setdefault("weights_only", False)
    return _orig_torch_load(*a, **k)
torch.load = _torch_load_compat

from multiome_label_transfer import (
    run_multi_model,
    run_gex_only_model,
    transfer_labels,
    finalize_output,
)


def label_transfer_from_preprocessed(gex_path, ref_path, input_id, atac_path=None,
                                     scanvi_model_tar=None, max_epochs=None, batch_size=None,
                                     output_max_probability=False):
    atac_present = bool(atac_path)
    print(f"Mode: {'multiome (GEX + ATAC)' if atac_present else 'GEX-only (no ATAC)'}", flush=True)

    # ── 1. Load preprocessed inputs (no re-preprocessing) ────────────────────
    print("Loading preprocessed inputs (no re-preprocessing)...", flush=True)
    gex = sc.read_h5ad(gex_path)
    ref = sc.read_h5ad(ref_path)
    print(f"  GEX:       {gex.shape}", flush=True)
    print(f"  Reference: {ref.shape}", flush=True)
    atac_activity = None
    if atac_present:
        atac_activity = sc.read_h5ad(atac_path)
        print(f"  ATAC activity: {atac_activity.shape}", flush=True)

    timing = {}
    start = time.time()

    # ── 2. Obtain the SCANVI model: load a supplied one, else train ──────────
    if scanvi_model_tar:
        import scvi
        print(f"Loading supplied SCANVI model (skip training): {scanvi_model_tar}", flush=True)
        extract_dir = "scanvi_model_in"
        with tarfile.open(scanvi_model_tar) as tf:
            tf.extractall(extract_dir)
        model_dir = os.path.dirname(glob.glob(f"{extract_dir}/**/model.pt", recursive=True)[0])
        # Rebuild the concat exactly as training does (so genes/labels line up), align to the model's
        # saved gene set, and load. No training is performed.
        if atac_present:
            adatas, keys = [gex, atac_activity, ref], ["rna_unannotated", "atac_unannotated", "rna_annotated"]
        else:
            adatas, keys = [gex, ref], ["rna_unannotated", "rna_annotated"]
        data = ad.concat(adatas, join="inner", label="modality", keys=keys, index_unique="_")
        sc.pp.filter_genes(data, min_cells=5)
        data.obs["celltype_scanvi"] = "Unknown"
        scvi.model.SCANVI.prepare_query_anndata(data, model_dir)   # subset/pad to the model's var_names
        lvae = scvi.model.SCANVI.load(model_dir, adata=data)
        timing["Model Load"] = time.time() - start
        print(f"  Model loaded in {timing['Model Load']:.1f}s", flush=True)
    else:
        print("Training SCVI and SCANVI models...", flush=True)
        train_kwargs = {}
        if max_epochs is not None:
            train_kwargs["max_epochs"] = max_epochs
        if batch_size is not None:
            train_kwargs["batch_size"] = batch_size
        print(f"  training kwargs: {train_kwargs}", flush=True)
        if atac_present:
            data, vae, lvae = run_multi_model(gex, atac_activity, ref, **train_kwargs)
        else:
            data, vae, lvae = run_gex_only_model(gex, ref, **train_kwargs)
        timing["Model Training"] = time.time() - start
        print(f"  Model training complete in {timing['Model Training']:.1f}s", flush=True)

    # ── 3. Transfer labels using the trained/loaded SCANVI model ─────────────
    print("Transferring labels with SCANVI...", flush=True)
    start = time.time()
    _plot, data = transfer_labels(data, lvae)
    timing["Label Transfer"] = time.time() - start

    # ── 3b. Optional per-cell label confidence (max posterior probability) ───
    if output_max_probability:
        print("Computing per-cell max label probability (SCANVI predict soft=True)...", flush=True)
        soft_probs = lvae.predict(data, soft=True)
        data.obs["max_probability"] = soft_probs.max(axis=1).astype("float32")

    # ── 4. Propagate predicted labels back to the original matrices ──────────
    print("Propagating predicted labels to GEX matrix...", flush=True)
    gex.obs["final_annotation"] = data.obs.loc[gex.obs_names + "_rna_unannotated"]["C_scANVI"].to_numpy()
    if output_max_probability:
        gex.obs["max_probability"] = data.obs.loc[gex.obs_names + "_rna_unannotated"]["max_probability"].to_numpy()
    if atac_present:
        print("Propagating predicted labels to ATAC matrix...", flush=True)
        atac_activity.obs["final_annotation"] = data.obs.loc[
            atac_activity.obs_names + "_atac_unannotated"]["C_scANVI"].to_numpy()
        if output_max_probability:
            atac_activity.obs["max_probability"] = data.obs.loc[
                atac_activity.obs_names + "_atac_unannotated"]["max_probability"].to_numpy()

    # ── 5. Write annotated GEX (and ATAC) matrices ───────────────────────────
    print(f"Writing annotated matrices (input_id={input_id})...", flush=True)
    gex.write(f"{input_id}_gex_annotated_matrix.h5ad")
    print(f"  {input_id}_gex_annotated_matrix.h5ad:  {gex.shape}", flush=True)
    if atac_present:
        atac_activity.write(f"{input_id}_atac_annotated_matrix.h5ad")
        print(f"  {input_id}_atac_annotated_matrix.h5ad: {atac_activity.shape}", flush=True)

    # ── 6. Finalize SCANVI predictions (metadata + raw counts for SCP) ───────
    print("Finalizing SCANVI predictions...", flush=True)
    final_data = finalize_output(data)
    final_data.write(f"{input_id}_SCANVI_predictions.h5ad")
    print(f"  {input_id}_SCANVI_predictions.h5ad:    {final_data.shape}", flush=True)

    # ── 7. Emit the minimal SCANVI model (no adata), tar'd for reuse ─────────
    print("Saving minimal SCANVI model (no adata)...", flush=True)
    lvae.save("scanvi_model_out", save_anndata=False, overwrite=True)
    with tarfile.open(f"{input_id}_scanvi_model.tar.gz", "w:gz") as tf:
        tf.add("scanvi_model_out")

    print("Timing:", {k: f"{v:.1f}s" for k, v in timing.items()}, flush=True)
    print("label_transfer_from_preprocessed complete.", flush=True)


def _parse_args():
    p = argparse.ArgumentParser(description="scANVI label transfer from preprocessed inputs.")
    p.add_argument("--gex", required=True, help="Preprocessed GEX h5ad.")
    p.add_argument("--ref", required=True, help="Preprocessed reference h5ad.")
    p.add_argument("--input-id", required=True, help="Prefix for all output filenames.")
    p.add_argument("--atac", default=None, help="Preprocessed ATAC activity h5ad (multiome mode).")
    p.add_argument("--scanvi-model", default=None, help="Saved SCANVI model .tar.gz to load (skip training).")
    p.add_argument("--max-epochs", type=int, default=None, help="SCVI/SCANVI training epochs (container default 500).")
    p.add_argument("--batch-size", type=int, default=None, help="SCVI/SCANVI minibatch size (container default 128).")
    p.add_argument("--output-max-probability", action="store_true",
                   help="Also write a max_probability obs column (the assigned label's confidence).")
    return p.parse_args()


if __name__ == "__main__":
    a = _parse_args()
    label_transfer_from_preprocessed(
        gex_path=a.gex, ref_path=a.ref, input_id=a.input_id, atac_path=a.atac,
        scanvi_model_tar=a.scanvi_model, max_epochs=a.max_epochs, batch_size=a.batch_size,
        output_max_probability=a.output_max_probability,
    )
