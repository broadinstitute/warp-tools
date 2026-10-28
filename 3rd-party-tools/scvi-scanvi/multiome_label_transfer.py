"""SCVI/SCANVI training and label transfer, imported by label_transfer_from_preprocessed.py.

  run_multi_model(gex, atac_activity, ref)   train on GEX + ATAC gene activity + reference
  run_gex_only_model(gex, ref)               train on GEX + reference
  transfer_labels(data, lvae)                predict labels (C_scANVI), latent (X_scANVI), UMAP
  finalize_output(data)                      rename final_annotation -> celltype, counts into .raw

Both training entry points share _concat_and_train_models, so the two modes cannot drift apart.
Inputs are preprocessed by preprocess_query.py (the query) and the WDL's PreprocessReference task
(the reference). Nothing here writes files; the caller writes the outputs.
"""
import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import scanpy as sc
import scvi


def _concat_and_train_models(adatas, keys, max_epochs=500, batch_size=128):
    """Concatenate the datasets and train SCVI, then SCANVI from it. Returns (data, vae, lvae).

    adatas, keys: datasets and their modality keys, in matching order. Each dataset needs a 'batch'
    column; the reference also needs 'final_annotation'. keys go into obs['modality'] and onto
    obs_names (index_unique='_'), e.g. '<barcode>_rna_unannotated', so labels can be copied back.
    max_epochs: cap for both models (SCVI also stops early). batch_size: SGD minibatch size; SCANVI
    activation memory scales with batch_size x number of labels.
    """
    data = ad.concat(adatas, join="inner", label="modality", keys=keys, index_unique="_")
    sc.pp.filter_genes(data, min_cells=5)
    sc.pp.highly_variable_genes(data, n_top_genes=5000, flavor="seurat_v3", batch_key="batch", subset=True)

    scvi.model.SCVI.setup_anndata(data, batch_key="batch")
    vae = scvi.model.SCVI(data, n_layers=2, n_latent=30, gene_likelihood="nb", dispersion="gene-batch")
    vae.train(max_epochs=max_epochs, batch_size=batch_size, early_stopping=True,
              datasplitter_kwargs={"drop_last": True})

    # Reference cells keep their labels; every query cell is 'Unknown' (unlabelled) for SCANVI.
    data.obs["celltype_scanvi"] = "Unknown"
    ref_idx = data.obs["modality"] == "rna_annotated"
    data.obs.loc[ref_idx, "celltype_scanvi"] = data.obs.loc[ref_idx, "final_annotation"]
    lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=data, labels_key="celltype_scanvi",
                                            unlabeled_category="Unknown")
    lvae.train(max_epochs=max_epochs, batch_size=batch_size, n_samples_per_label=100,
               datasplitter_kwargs={"drop_last": True})
    return data, vae, lvae


def run_multi_model(gex, atac_activity_matrix, ref, max_epochs=500, batch_size=128):
    """Train on GEX + ATAC gene activity + reference (multiome mode)."""
    return _concat_and_train_models([gex, atac_activity_matrix, ref],
                                    keys=["rna_unannotated", "atac_unannotated", "rna_annotated"],
                                    max_epochs=max_epochs, batch_size=batch_size)


def run_gex_only_model(gex, ref, max_epochs=500, batch_size=128):
    """Train on GEX + reference (GEX-only mode); the GEX key is the same as in multiome mode."""
    return _concat_and_train_models([gex, ref], keys=["rna_unannotated", "rna_annotated"],
                                    max_epochs=max_epochs, batch_size=batch_size)


def transfer_labels(data, lvae):
    """Predict a label for every cell (obs['C_scANVI']), store the latent (obsm['X_scANVI']), and
    fit neighbors and a UMAP on that latent. Returns data."""
    data.obs["C_scANVI"] = lvae.predict(data)
    data.obsm["X_scANVI"] = lvae.get_latent_representation(data)
    sc.pp.neighbors(data, use_rep="X_scANVI")
    sc.tl.umap(data)
    return data


def finalize_output(scanvi_output):
    """Rename final_annotation to celltype and copy the counts into .raw. Adds no sample metadata."""
    scanvi_output.obs.rename(columns={"final_annotation": "celltype"}, inplace=True)
    scanvi_output.raw = ad.AnnData(scanvi_output.X, var=scanvi_output.var, obs=scanvi_output.obs)
    return scanvi_output
