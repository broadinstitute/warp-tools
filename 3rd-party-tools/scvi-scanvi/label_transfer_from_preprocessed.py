"""scANVI label transfer for one sample, after preprocessing.

Labels the sample at every label level of the model (model_archive.py): each level has its own
SCANVI model, and each gives the sample a label column, a latent and a UMAP. The models come from
one of three places:
  - no model supplied: train SCVI on the reference plus this sample, then one SCANVI per level of
    the reference (needs --ref);
  - a model this pipeline trained: it already knows the sample's batch label, so each level loads
    and predicts on the sample alone;
  - a model trained on the reference alone (warp-scanvi-training): the sample's batch is new to
    it, so the sample is mapped onto each level with scArches (load_query_data plus a short
    training on the sample alone).
With a supplied model the reference is not needed. A supplied model whose genes are Ensembl IDs is
matched through the sample's Ensembl IDs, looked up from its Optimus ensembl_ids column.

The annotated matrices are written SCP-ready (see scp_output.py) and record the run provenance in
uns: the sample's id and input file, the model's name (the reference it was trained on), its
levels and the image's library versions. No sample metadata is added. scvi-tools picks the
accelerator itself, so the same command runs on GPU and CPU machines.

  python3 label_transfer_from_preprocessed.py --gex g.h5ad --input-id ID --gex-out g_annotated.h5ad \
      (--scanvi-model m.tar.gz | --ref r.h5ad [--model-out m.tar.gz]) [--atac a.h5ad --atac-out ...] \
      [--max-epochs N] [--query-max-epochs N] [--batch-size N] [--seed N] [--output-max-probability]
"""
import argparse
import importlib.metadata
import os
import re
import time

import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import torch

from model_archive import read_archive, write_archive
from multiome_label_transfer import train_scanvi, train_scvi
from scp_output import make_scp_ready, write_h5ad

QUERY_MAX_EPOCHS = 100   # scvi-tools' scArches tutorial value for training on the query
ENSEMBL_ID = re.compile(r"ENS[A-Z]*G\d+")


def match_gene_ids(data, gex, model_dir):
    """Name the sample's genes (data, the sample's modelling copy) the way the model does: for a
    model whose genes are Ensembl IDs, each gene symbol is looked up in the sample's Ensembl ID
    column. The annotated outputs keep the sample's own names. Prints the genes in common."""
    model_genes = list(torch.load(os.path.join(model_dir, "model.pt"), map_location="cpu",
                                  weights_only=False)["var_names"])
    if sum(bool(ENSEMBL_ID.fullmatch(g)) for g in model_genes[:100]) > 50:
        column = next((c for c in ("ensembl_ids", "gene_ids", "ensembl_id") if c in gex.var), None)
        if column is None:
            raise SystemExit("ERROR: the model names genes by Ensembl ID, but the sample has no Ensembl ID "
                             "column (Optimus writes var['ensembl_ids']).")
        lookup = dict(zip(gex.var_names, gex.var[column].astype(str).str.replace(r"\.\d+$", "", regex=True)))
        data.var_names = [lookup.get(g, g) for g in data.var_names]
        print(f"  The model names genes by Ensembl ID; looked them up in var['{column}'].", flush=True)
    found = len(set(model_genes) & set(data.var_names))
    print(f"  The sample has {found:,} of the model's {len(model_genes):,} genes.", flush=True)
    return data


def load_supplied_model(data, model_dir, query_max_epochs=None, batch_size=None):
    """Fit one level's saved SCANVI model in model_dir to this sample alone. Returns (query, model),
    where query is data with its genes matched to the model's."""
    query = scvi.model.SCANVI.prepare_query_anndata(data, model_dir, inplace=False)
    registry = scvi.model.SCANVI.load_registry(model_dir)
    known = registry["field_registries"]["batch"]["state_registry"]["categorical_mapping"]
    new = sorted(set(query.obs["batch"].astype(str)) - set(map(str, known)))
    if not new:
        print("  The model knows this sample's batch label; loading it as is.", flush=True)
        return query, scvi.model.SCANVI.load(model_dir, adata=query)

    epochs = query_max_epochs or QUERY_MAX_EPOCHS
    print(f"  Batch {new} is new to the model; mapping the sample onto it (scArches, "
          f"{epochs} epochs on the sample).", flush=True)
    model = scvi.model.SCANVI.load_query_data(query, model_dir)
    model.train(max_epochs=epochs, batch_size=batch_size or 128,
                plan_kwargs={"weight_decay": 0.0}, datasplitter_kwargs={"drop_last": True})
    return query, model


def _predict(lvae, adata, with_probability):
    """(labels, latent, max posterior or None) for every cell of adata."""
    labels = np.asarray(lvae.predict(adata)).astype(str)
    latent = lvae.get_latent_representation(adata)
    probability = None
    if with_probability:
        probability = lvae.predict(adata, soft=True).max(axis=1).to_numpy(dtype="float32")
    return labels, latent, probability


def _annotate(matrix, obs_names, suffix, results, provenance):
    """Copy this matrix's rows of each level's labels and latent, then make it SCP-ready."""
    rows = obs_names.get_indexer(matrix.obs_names + suffix)
    if (rows < 0).any():
        raise ValueError(f"{(rows < 0).sum()} cells of the matrix are missing from the labelled data.")
    for level, (labels, _, probability) in results.items():
        if level in matrix.obs:
            print(f"  The input's obs['{level}'] is replaced by the predicted labels.", flush=True)
        matrix.obs[level] = labels[rows]
        if probability is not None:
            matrix.obs[f"{level}_max_probability"] = probability[rows]
    latents = {level: latent[rows] for level, (_, latent, _) in results.items()}
    return make_scp_ready(matrix, latents, provenance)


def label_transfer(gex, input_id, gex_out, ref=None, atac_activity=None, atac_out=None,
                   scanvi_model_tar=None, model_out=None, max_epochs=None, query_max_epochs=None,
                   batch_size=None, output_max_probability=False, seed=None, input_file=""):
    """Label one preprocessed sample (in memory) at every level and write its outputs."""
    if not scanvi_model_tar and ref is None:
        raise SystemExit("ERROR: give a reference (to train) or a SCANVI model (to predict).")
    # Unset (default): training is nondeterministic. Tests pin a seed to make runs reproducible.
    if seed is not None:
        scvi.settings.seed = seed

    results = {}   # level -> (labels, latent, max posterior), one row per cell of data
    start = time.time()
    if scanvi_model_tar:
        print(f"Using the supplied model {scanvi_model_tar}", flush=True)
        reference, level_dirs = read_archive(scanvi_model_tar)
        adatas, keys = [gex], ["rna_unannotated"]
        if atac_activity is not None:
            adatas, keys = [gex, atac_activity], ["rna_unannotated", "atac_unannotated"]
        data = ad.concat(adatas, join="inner", label="modality", keys=keys, index_unique="_")
        data.obs["celltype_scanvi"] = "Unknown"
        data = match_gene_ids(data, gex, next(iter(level_dirs.values())))
        for level, model_dir in level_dirs.items():
            print(f"Level '{level}':", flush=True)
            query, lvae = load_supplied_model(data, model_dir, query_max_epochs, batch_size)
            results[level] = _predict(lvae, query, output_max_probability)
    else:
        reference = ref.uns.get("scanvi_reference")   # set by the WDL's PreprocessReference
        levels = [str(level) for level in ref.uns["scanvi_levels"]]
        train_kwargs = {k: v for k, v in [("max_epochs", max_epochs), ("batch_size", batch_size)] if v is not None}
        print("Training SCVI on the reference plus this sample...", flush=True)
        data, vae = train_scvi(gex, atac_activity, ref, **train_kwargs)
        if model_out:
            vae.save("models_out/scvi", save_anndata=False, overwrite=True)
        for level in levels:
            print(f"Training SCANVI on level '{level}'...", flush=True)
            lvae = train_scanvi(data, vae, ref, level, **train_kwargs)
            results[level] = _predict(lvae, data, output_max_probability)
            if model_out:
                lvae.save(f"models_out/scanvi/{level}", save_anndata=False, overwrite=True)
        if model_out:
            write_archive(model_out, str(reference), "models_out/scvi",
                          {level: f"models_out/scanvi/{level}" for level in levels})
    print(f"  Labelled at levels {list(results)} in {time.time() - start:.1f}s", flush=True)

    # scanvi_model names the model by the reference it was trained on (reference.txt); a model
    # without one falls back to its archive's file name.
    provenance = {
        "scanvi_input_id": input_id,
        "scanvi_input_file": input_file,
        "scanvi_model": str(reference) if reference else os.path.basename(scanvi_model_tar or model_out or ""),
        "scanvi_levels": list(results),
        "scanvi_versions": {p: importlib.metadata.version(p) for p in
                            ("scvi-tools", "torch", "scanpy", "anndata", "snapatac2", "numpy", "scipy")},
    }
    print(f"Writing outputs (input_id={input_id})...", flush=True)
    gex = _annotate(gex, data.obs_names, "_rna_unannotated", results, provenance)
    write_h5ad(gex, gex_out)
    if atac_activity is not None:
        atac_activity = _annotate(atac_activity, data.obs_names, "_atac_unannotated", results, provenance)
        write_h5ad(atac_activity, atac_out)
    print(f"Done: {input_id} ({gex.n_obs} GEX cells) -> {os.path.basename(gex_out)}", flush=True)


def add_label_transfer_args(p):
    p.add_argument("--input-id", required=True, help="The sample's id, recorded in uns['scanvi_input_id'].")
    p.add_argument("--gex-out", required=True, help="Path of the annotated GEX h5ad to write.")
    p.add_argument("--atac-out", default=None, help="Path of the annotated ATAC h5ad to write (multiome).")
    p.add_argument("--ref", default=None, help="Preprocessed reference h5ad (needed to train).")
    p.add_argument("--scanvi-model", default=None, help="Model archive .tar.gz to use instead of training.")
    p.add_argument("--model-out", default=None, help="When training, also write the trained models here.")
    p.add_argument("--max-epochs", type=int, default=None, help="SCVI/SCANVI training epochs (container default 500).")
    p.add_argument("--query-max-epochs", type=int, default=None,
                   help=f"Epochs when mapping a sample onto a reference-only model (default {QUERY_MAX_EPOCHS}).")
    p.add_argument("--batch-size", type=int, default=None, help="SCVI/SCANVI minibatch size (container default 128).")
    p.add_argument("--output-max-probability", action="store_true",
                   help="Also write a <level>_max_probability obs column (each label's confidence).")
    p.add_argument("--seed", type=int, default=None,
                   help="Random seed for reproducible SCVI/SCANVI training (default: unseeded).")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="scANVI label transfer from preprocessed inputs.")
    p.add_argument("--gex", required=True, help="Preprocessed GEX h5ad.")
    p.add_argument("--atac", default=None, help="Preprocessed ATAC activity h5ad (multiome mode).")
    add_label_transfer_args(p)
    a = p.parse_args()
    label_transfer(sc.read_h5ad(a.gex), a.input_id, a.gex_out,
                   ref=sc.read_h5ad(a.ref) if a.ref and not a.scanvi_model else None,
                   atac_activity=sc.read_h5ad(a.atac) if a.atac else None, atac_out=a.atac_out,
                   scanvi_model_tar=a.scanvi_model, model_out=a.model_out, max_epochs=a.max_epochs,
                   query_max_epochs=a.query_max_epochs, batch_size=a.batch_size,
                   output_max_probability=a.output_max_probability, seed=a.seed,
                   input_file=os.path.basename(a.gex))
