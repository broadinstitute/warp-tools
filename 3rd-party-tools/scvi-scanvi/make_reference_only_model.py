"""Train a SCANVI model on a reference alone, the way warp-scanvi-training does.

Builds the scANVI pipeline's query-mapping test fixture, and is used by the image's CI test. Same
setup as the harness (rsc3/warp-scanvi-training, settings of 2026-10-29): every gene of the
reference; batch -> 'batch' and labels -> 'celltype_scanvi' as plain strings, 'Unknown' unlabeled;
scVI set up without labels, with the harness's SCVI_MODEL_KWARGS on 2 layers and a 30-dimensional
latent; then SCANVI from it. Saved without data and tarred, the form the pipeline's scanvi_model
input takes, with the reference's id (reference_id.py) in reference.txt. The scArches settings
converge slowly on a small reference: about 200 epochs for an 800-cell one.

  python3 make_reference_only_model.py --ref ref.h5ad --out model.tar.gz \
      [--batch-key donor_id] [--labels-key subclass] [--epochs 200] [--seed 0]
"""
import argparse
import tarfile

import anndata as ad
import scvi

from reference_id import reference_id, write_reference

# Keep identical to SCVI_MODEL_KWARGS in warp-scanvi-training's train_atlas.py: negative binomial;
# dispersion 'gene'; scvi-tools' scArches settings, so samples can be mapped onto the model.
SCVI_MODEL_KWARGS = {"gene_likelihood": "nb", "dispersion": "gene", "encode_covariates": True,
                     "use_layer_norm": "both", "use_batch_norm": "none", "dropout_rate": 0.2}


def train_reference_only_model(ref, batch_key, labels_key, epochs, seed, out_tar, reference=None):
    if ref.X is None:   # AIT references keep their counts in .raw
        ref = ad.AnnData(X=ref.raw.X, obs=ref.obs, var=ref.var)
    ref = ref.copy()
    ref.obs["batch"] = ref.obs[batch_key].astype(str)
    ref.obs["celltype_scanvi"] = ref.obs[labels_key].astype(str)
    scvi.settings.seed = seed
    scvi.model.SCVI.setup_anndata(ref, batch_key="batch")
    vae = scvi.model.SCVI(ref, n_layers=2, n_latent=30, **SCVI_MODEL_KWARGS)
    vae.train(max_epochs=epochs, datasplitter_kwargs={"drop_last": True})
    lvae = scvi.model.SCANVI.from_scvi_model(vae, labels_key="celltype_scanvi", unlabeled_category="Unknown")
    lvae.train(max_epochs=epochs, datasplitter_kwargs={"drop_last": True})
    lvae.save("reference_only_model", save_anndata=False, overwrite=True)
    if reference:
        write_reference("reference_only_model", reference)
    with tarfile.open(out_tar, "w:gz") as tf:
        tf.add("reference_only_model")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--batch-key", default="donor_id")
    p.add_argument("--labels-key", default="subclass")
    p.add_argument("--epochs", type=int, default=200)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    ref = ad.read_h5ad(a.ref)
    train_reference_only_model(ref, a.batch_key, a.labels_key, a.epochs, a.seed, a.out,
                               reference=reference_id(a.ref, ref.uns))
