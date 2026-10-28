"""Put an annotated matrix in the form Single Cell Portal (SCP) reads.

make_scp_ready(adata, latent, provenance) changes the matrix in place:
  .raw              the counts (float32), over every gene the matrix keeps
  .X                counts normalized to 10,000 per cell, then log1p (float32)
  obsm['X_scANVI']  the full scANVI latent for these cells
  obsm['X_umap']    a 2-column UMAP fitted on that latent, for these cells only
  uns               run provenance: scanvi_input_id, scanvi_input_file, scanvi_model
It adds no sample metadata. write_h5ad() writes with gzip compression.
"""
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

TARGET_SUM = 1e4


def make_scp_ready(adata, latent, provenance):
    counts = sp.csr_matrix(adata.layers.pop("counts", adata.X), dtype=np.float32)
    # .raw gets its own copy: normalize_total and log1p change .X in place.
    adata.raw = ad.AnnData(X=counts.copy(), obs=pd.DataFrame(index=adata.obs_names), var=adata.var.copy())
    adata.X = counts
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    # scanpy records base=None (natural log), which anndata before 0.11 cannot read back; an empty
    # entry means the same to scanpy and keeps the "already log-transformed" marker.
    adata.uns["log1p"] = {}
    adata.obsm["X_scANVI"] = np.asarray(latent, dtype=np.float32)
    sc.pp.neighbors(adata, use_rep="X_scANVI", key_added="scp")
    sc.tl.umap(adata, neighbors_key="scp", n_components=2)
    # The neighbor graph only served the UMAP; SCP does not read it.
    del adata.obsp["scp_connectivities"], adata.obsp["scp_distances"]
    del adata.uns["scp"], adata.uns["umap"]
    adata.uns.update(provenance)
    return adata


def write_h5ad(adata, path):
    adata.write_h5ad(path, compression="gzip")
