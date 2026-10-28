"""scANVI label transfer for one sample, after preprocessing.

Obtains a SCANVI model, transfers labels, optionally records the per-cell max posterior, and writes
the outputs. The model comes from one of three places:
  - no model supplied: train SCVI + SCANVI on the reference plus this sample (needs --ref);
  - a model this pipeline trained: it already knows the sample's batch label, so it loads and
    predicts on the sample alone;
  - a model trained on the reference alone (warp-scanvi-training): the sample's batch is new to
    it, so the sample is mapped onto it with scArches (load_query_data plus a short training on
    the sample alone).
With a supplied model the reference is not needed and the outputs hold the sample's cells only.

The annotated matrices are written SCP-ready (see scp_output.py) and every output records the run
provenance in uns: the sample's id and input file, and in scanvi_model the reference the model was
trained on (reference_id.py), carried in the model's reference.txt. No sample metadata is added.
scvi-tools picks the accelerator itself, so the same command runs on GPU and CPU machines.

Writes to the working directory:
  <input_id>_gex_annotated_matrix.h5ad
  <input_id>_atac_annotated_matrix.h5ad   (multiome only)
  <input_id>_SCANVI_predictions.h5ad
  <input_id>_scanvi_model.tar.gz          (minimal SCANVI model, no bundled adata)

  python3 label_transfer_from_preprocessed.py --gex g.h5ad --input-id ID \
      (--scanvi-model m.tar.gz | --ref r.h5ad) [--atac a.h5ad] [--max-epochs N] \
      [--query-max-epochs N] [--batch-size N] [--seed N] [--output-max-probability]
"""
import argparse
import glob
import importlib.metadata
import os
import tarfile
import time

import anndata as ad
import scanpy as sc
import scvi

from multiome_label_transfer import (
    run_multi_model,
    run_gex_only_model,
    transfer_labels,
    finalize_output,
)
from reference_id import read_reference, write_reference
from scp_output import make_scp_ready, write_h5ad

QUERY_MAX_EPOCHS = 100   # scvi-tools' scArches tutorial value for training on the query


def _extract_model(model_tar):
    with tarfile.open(model_tar) as tf:
        tf.extractall("scanvi_model_in", filter="data")
    matches = glob.glob("scanvi_model_in/**/model.pt", recursive=True)
    if len(matches) != 1:
        raise ValueError(f"Expected exactly 1 model.pt in archive, found {len(matches)}.")
    return os.path.dirname(matches[0])


def load_supplied_model(gex, atac_activity, model_dir, query_max_epochs=None, batch_size=None):
    """Fit the saved SCANVI model in model_dir to this sample alone. Returns (data, model)."""
    adatas, keys = [gex], ["rna_unannotated"]
    if atac_activity is not None:
        adatas, keys = [gex, atac_activity], ["rna_unannotated", "atac_unannotated"]
    data = ad.concat(adatas, join="inner", label="modality", keys=keys, index_unique="_")
    data.obs["celltype_scanvi"] = "Unknown"
    data = scvi.model.SCANVI.prepare_query_anndata(data, model_dir, inplace=False)

    registry = scvi.model.SCANVI.load_registry(model_dir)
    known = registry["field_registries"]["batch"]["state_registry"]["categorical_mapping"]
    new = sorted(set(data.obs["batch"].astype(str)) - set(map(str, known)))
    if not new:
        print("  The model knows this sample's batch label; loading it as is.", flush=True)
        return data, scvi.model.SCANVI.load(model_dir, adata=data)

    epochs = query_max_epochs or QUERY_MAX_EPOCHS
    print(f"  Batch {new} is new to the model; mapping the sample onto it (scArches, "
          f"{epochs} epochs on the sample).", flush=True)
    model = scvi.model.SCANVI.load_query_data(data, model_dir)
    model.train(max_epochs=epochs, batch_size=batch_size or 128,
                plan_kwargs={"weight_decay": 0.0}, datasplitter_kwargs={"drop_last": True})
    return data, model


def _annotate(matrix, data, suffix, provenance, with_probability):
    """Copy this matrix's labels and latent rows from the joint object, then make it SCP-ready."""
    rows = data.obs_names.get_indexer(matrix.obs_names + suffix)
    if (rows < 0).any():
        raise ValueError(f"{(rows < 0).sum()} cells of the matrix are missing from the labelled data.")
    matrix.obs["final_annotation"] = data.obs["C_scANVI"].to_numpy()[rows]
    if with_probability:
        matrix.obs["max_probability"] = data.obs["max_probability"].to_numpy()[rows]
    return make_scp_ready(matrix, data.obsm["X_scANVI"][rows], provenance)


def label_transfer(gex, input_id, ref=None, atac_activity=None, scanvi_model_tar=None,
                   max_epochs=None, query_max_epochs=None, batch_size=None,
                   output_max_probability=False, seed=None, input_file=""):
    """Label one preprocessed sample (in memory) and write its outputs."""
    if not scanvi_model_tar and ref is None:
        raise SystemExit("ERROR: give a reference (to train) or a SCANVI model (to predict).")
    # Unset (default): training is nondeterministic. Tests pin a seed to make runs reproducible.
    if seed is not None:
        scvi.settings.seed = seed

    start = time.time()
    if scanvi_model_tar:
        print(f"Using the supplied SCANVI model {scanvi_model_tar}", flush=True)
        model_dir = _extract_model(scanvi_model_tar)
        reference = read_reference(model_dir)
        data, lvae = load_supplied_model(gex, atac_activity, model_dir, query_max_epochs, batch_size)
    else:
        reference = ref.uns.get("scanvi_reference")   # set by the WDL's PreprocessReference
        print("Training SCVI and SCANVI on the reference plus this sample...", flush=True)
        train_kwargs = {k: v for k, v in [("max_epochs", max_epochs), ("batch_size", batch_size)] if v is not None}
        if atac_activity is not None:
            data, _vae, lvae = run_multi_model(gex, atac_activity, ref, **train_kwargs)
        else:
            data, _vae, lvae = run_gex_only_model(gex, ref, **train_kwargs)
    print(f"  Model ready in {time.time() - start:.1f}s", flush=True)

    print("Transferring labels with SCANVI...", flush=True)
    data = transfer_labels(data, lvae)
    if output_max_probability:
        data.obs["max_probability"] = lvae.predict(data, soft=True).max(axis=1).astype("float32")

    # The model is named by the reference it was trained on; a model without reference.txt falls
    # back to its archive's file name. scanvi_versions records the image's key library versions.
    provenance = {
        "scanvi_input_id": input_id,
        "scanvi_input_file": input_file,
        "scanvi_model": str(reference) if reference else (
            os.path.basename(scanvi_model_tar) if scanvi_model_tar else f"{input_id}_scanvi_model.tar.gz"),
        "scanvi_versions": {p: importlib.metadata.version(p) for p in
                            ("scvi-tools", "torch", "scanpy", "anndata", "snapatac2", "numpy", "scipy")},
    }
    print(f"Writing outputs (input_id={input_id})...", flush=True)
    gex = _annotate(gex, data, "_rna_unannotated", provenance, output_max_probability)
    write_h5ad(gex, f"{input_id}_gex_annotated_matrix.h5ad")
    if atac_activity is not None:
        atac_activity = _annotate(atac_activity, data, "_atac_unannotated", provenance, output_max_probability)
        write_h5ad(atac_activity, f"{input_id}_atac_annotated_matrix.h5ad")

    data.obs["final_annotation"] = data.obs["C_scANVI"]
    predictions = finalize_output(data)
    predictions.uns.update(provenance)
    write_h5ad(predictions, f"{input_id}_SCANVI_predictions.h5ad")

    lvae.save("scanvi_model_out", save_anndata=False, overwrite=True)
    if reference:
        write_reference("scanvi_model_out", str(reference))
    with tarfile.open(f"{input_id}_scanvi_model.tar.gz", "w:gz") as tf:
        tf.add("scanvi_model_out")
    print(f"Done: {input_id} ({gex.n_obs} GEX cells, predictions {predictions.shape}).", flush=True)


def label_transfer_from_preprocessed(gex_path, ref_path, input_id, atac_path=None,
                                     scanvi_model_tar=None, max_epochs=None, batch_size=None,
                                     output_max_probability=False, seed=None,
                                     query_max_epochs=None, input_file=None):
    """Path-based wrapper around label_transfer(); ref_path may be None when a model is supplied."""
    gex = sc.read_h5ad(gex_path)
    ref = sc.read_h5ad(ref_path) if ref_path and not scanvi_model_tar else None
    atac_activity = sc.read_h5ad(atac_path) if atac_path else None
    label_transfer(gex, input_id, ref=ref, atac_activity=atac_activity,
                   scanvi_model_tar=scanvi_model_tar, max_epochs=max_epochs,
                   query_max_epochs=query_max_epochs, batch_size=batch_size,
                   output_max_probability=output_max_probability, seed=seed,
                   input_file=input_file or os.path.basename(gex_path))


def add_label_transfer_args(p):
    p.add_argument("--input-id", required=True, help="Prefix for all output filenames.")
    p.add_argument("--ref", default=None, help="Preprocessed reference h5ad (needed to train).")
    p.add_argument("--scanvi-model", default=None, help="Saved SCANVI model .tar.gz to use instead of training.")
    p.add_argument("--max-epochs", type=int, default=None, help="SCVI/SCANVI training epochs (container default 500).")
    p.add_argument("--query-max-epochs", type=int, default=None,
                   help=f"Epochs when mapping a sample onto a reference-only model (default {QUERY_MAX_EPOCHS}).")
    p.add_argument("--batch-size", type=int, default=None, help="SCVI/SCANVI minibatch size (container default 128).")
    p.add_argument("--output-max-probability", action="store_true",
                   help="Also write a max_probability obs column (the assigned label's confidence).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducible SCVI/SCANVI training (default: unseeded).")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="scANVI label transfer from preprocessed inputs.")
    p.add_argument("--gex", required=True, help="Preprocessed GEX h5ad.")
    p.add_argument("--atac", default=None, help="Preprocessed ATAC activity h5ad (multiome mode).")
    add_label_transfer_args(p)
    a = p.parse_args()
    label_transfer_from_preprocessed(
        gex_path=a.gex, ref_path=a.ref, input_id=a.input_id, atac_path=a.atac,
        scanvi_model_tar=a.scanvi_model, max_epochs=a.max_epochs, batch_size=a.batch_size,
        output_max_probability=a.output_max_probability, seed=a.seed,
        query_max_epochs=a.query_max_epochs,
    )
