"""SCVI/SCANVI training on a reference plus one sample, imported by label_transfer_from_preprocessed.py.

  train_scvi(gex, atac_activity, ref)       concatenate and train SCVI; returns (data, vae)
  train_scanvi(data, vae, ref, level)       train SCANVI from vae on one label level of the reference

scVI is trained once; each label level then gets its own SCANVI model from it. Inputs are
preprocessed by preprocess_query.py (the query) and the WDL's PreprocessReference task (the
reference, which lists its levels in uns['scanvi_levels']). Nothing here writes files.
"""
import warnings
warnings.filterwarnings("ignore")

import anndata as ad
import scanpy as sc
import scvi

REF_KEY = "rna_annotated"


def train_scvi(gex, atac_activity, ref, max_epochs=500, batch_size=128):
    """Concatenate the sample (GEX, plus ATAC gene activity in multiome mode) with the reference and
    train SCVI. Returns (data, vae).

    Each dataset needs a 'batch' column. The modality keys go into obs['modality'] and onto obs_names
    (index_unique='_'), e.g. '<barcode>_rna_unannotated', so labels can be copied back. max_epochs:
    cap for SCVI (which also stops early) and each SCANVI. batch_size: SGD minibatch size; SCANVI
    activation memory scales with batch_size x number of labels.
    """
    adatas, keys = [gex, ref], ["rna_unannotated", REF_KEY]
    if atac_activity is not None:
        adatas, keys = [gex, atac_activity, ref], ["rna_unannotated", "atac_unannotated", REF_KEY]
    data = ad.concat(adatas, join="inner", label="modality", keys=keys, index_unique="_")
    sc.pp.filter_genes(data, min_cells=5)
    sc.pp.highly_variable_genes(data, n_top_genes=5000, flavor="seurat_v3", batch_key="batch", subset=True)

    scvi.model.SCVI.setup_anndata(data, batch_key="batch")
    vae = scvi.model.SCVI(data, n_layers=2, n_latent=30, gene_likelihood="nb", dispersion="gene-batch")
    vae.train(max_epochs=max_epochs, batch_size=batch_size, early_stopping=True,
              datasplitter_kwargs={"drop_last": True})
    return data, vae


def train_scanvi(data, vae, ref, level, max_epochs=500, batch_size=128):
    """Train SCANVI from vae with the reference's obs[level] as labels. Every sample cell is
    'Unknown' (unlabelled). Use the model before training the next level: each level's labels
    replace the last in data.obs['celltype_scanvi']."""
    labels = ref.obs[level].astype(str)
    labels.index = ref.obs_names + "_" + REF_KEY
    data.obs["celltype_scanvi"] = labels.reindex(data.obs_names).fillna("Unknown").to_numpy()
    lvae = scvi.model.SCANVI.from_scvi_model(vae, adata=data, labels_key="celltype_scanvi",
                                            unlabeled_category="Unknown")
    lvae.train(max_epochs=max_epochs, batch_size=batch_size, n_samples_per_label=100,
               datasplitter_kwargs={"drop_last": True})
    return lvae
