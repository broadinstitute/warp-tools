"""Label one query sample end to end: preprocess, transfer labels, write outputs.

scANVI.wdl's label-transfer task runs this once per sample, each in its own process and working
directory, so memory is freed between samples. Outputs are those of
label_transfer_from_preprocessed.label_transfer(), named by --input-id.

  python3 scanvi_sample.py --input-id ID --gex g.h5ad [--atac a.h5ad] [--genome hg38] \
      (--scanvi-model m.tar.gz | --ref preprocessed_ref.h5ad) [--max-epochs N] \
      [--query-max-epochs N] [--batch-size N] [--seed N] [--output-max-probability]
"""
import argparse
import os

import scanpy as sc

from label_transfer_from_preprocessed import add_label_transfer_args, label_transfer
from preprocess_query import preprocess_query

if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gex", required=True, help="The sample's GEX h5ad (Optimus output).")
    p.add_argument("--atac", default=None, help="The sample's ATAC cell-by-bin h5ad (multiome mode).")
    p.add_argument("--genome", default="hg38", help="Genome for the ATAC gene-activity conversion.")
    add_label_transfer_args(p)
    a = p.parse_args()

    print(f"=== {a.input_id}: preprocessing ===", flush=True)
    gex, atac_activity = preprocess_query(a.gex, a.atac, a.genome)
    ref = sc.read_h5ad(a.ref) if a.ref and not a.scanvi_model else None
    print(f"=== {a.input_id}: label transfer ===", flush=True)
    label_transfer(gex, a.input_id, ref=ref, atac_activity=atac_activity,
                   scanvi_model_tar=a.scanvi_model, max_epochs=a.max_epochs,
                   query_max_epochs=a.query_max_epochs, batch_size=a.batch_size,
                   output_max_probability=a.output_max_probability, seed=a.seed,
                   input_file=os.path.basename(a.gex))
