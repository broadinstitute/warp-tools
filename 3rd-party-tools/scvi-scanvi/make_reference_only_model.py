"""Train scANVI models on a reference alone, the way warp-scanvi-training does.

Builds the scANVI pipeline's query-mapping test fixture, and is used by the image's CI test. Same
setup as the harness (rsc3/warp-scanvi-training, settings of 2026-10-29): every gene of the
reference; batch -> 'batch' and each level's labels -> 'celltype_scanvi' as plain strings,
'Unknown' unlabeled; scVI set up without labels, with the harness's SCVI_MODEL_KWARGS on 2 layers
and a 30-dimensional latent; then one SCANVI from it per label level. Saved without data into the
model archive (model_archive.py), the form the pipeline's scanvi_model input takes, named by the
reference's id. The scArches settings converge slowly on a small reference: about 200 epochs for
an 800-cell one.

  python3 make_reference_only_model.py --ref ref.h5ad --out model.tar.gz \
      [--batch-key donor_id] [--levels class subclass] [--epochs 200] [--seed 0]
"""
import argparse

import anndata as ad
import scvi

from model_archive import reference_id, write_archive

# Keep identical to SCVI_MODEL_KWARGS in warp-scanvi-training's train_atlas.py: negative binomial;
# dispersion 'gene'; scvi-tools' scArches settings, so samples can be mapped onto the model.
SCVI_MODEL_KWARGS = {"gene_likelihood": "nb", "dispersion": "gene", "encode_covariates": True,
                     "use_layer_norm": "both", "use_batch_norm": "none", "dropout_rate": 0.2}


def train_reference_only_model(ref, batch_key, levels, epochs, seed, out_tar, reference):
    if ref.X is None:   # AIT references keep their counts in .raw
        ref = ad.AnnData(X=ref.raw.X, obs=ref.obs, var=ref.var)
    ref = ref.copy()
    ref.obs["batch"] = ref.obs[batch_key].astype(str)
    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(ref, batch_key="batch")
    vae = scvi.model.SCVI(ref, n_layers=2, n_latent=30, **SCVI_MODEL_KWARGS)
    vae.train(max_epochs=epochs, datasplitter_kwargs={"drop_last": True})
    vae.save("reference_only_model/scvi", save_anndata=False, overwrite=True)
    for level in levels:
        ref.obs["celltype_scanvi"] = ref.obs[level].astype(str)
        lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="celltype_scanvi", unlabeled_category="Unknown")
        lvae.train(max_epochs=epochs, datasplitter_kwargs={"drop_last": True})
        lvae.save(f"reference_only_model/scanvi/{level}", save_anndata=False, overwrite=True)
    write_archive(out_tar, reference, "reference_only_model/scvi",
                  {level: f"reference_only_model/scanvi/{level}" for level in levels})


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-key", default="donor_id")
    p.add_argument("--levels", nargs="+", default=["subclass"], help="Label columns, one SCANVI each.")
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    ref = ad.read_h5ad(a.ref)
    train_reference_only_model(ref, a.batch_key, a.levels, a.epochs, a.seed, a.out,
                               reference=reference_id(a.ref, ref.uns))
