# scvi-scanvi

## Quick reference

```bash
docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:635d4391d50cba9bd58f1fc41b10d8e1c61285a73bde75371815ce9a0db3430c
```

- **What is this image:** a GPU-enabled Python image (Debian `python:3.12-slim` base) for single-cell cell-type **label transfer** with deep generative models, bundling the `multiome_label_transfer.py` workflow.
- **What is scvi-tools:** [scvi-tools](https://docs.scvi-tools.org/) provides SCVI (unsupervised variational inference) and SCANVI (semi-supervised annotation) models for single-cell omics.

## Synopsis

This image trains SCVI/SCANVI models to transfer cell-type labels from an annotated reference onto unlabeled query cells. It supports two modes: **multiome** (gene expression + ATAC gene-activity + reference) and **GEX-only** (gene expression + reference, no ATAC). It is the execution environment for the Broad WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi); the workflow follows the [SnapATAC2 annotation tutorial](https://kzhang.org/SnapATAC2/tutorials/annotation.html) and uses the annotated PBMC reference from the scArches tutorial.

## Image contents

- Base image: `python:3.12-slim-trixie` (`--platform=linux/amd64`)
- `scvi-tools` 1.2 · `snapatac2` 2.7 · `scanpy` · `anndata` · `numpy` · `scikit-misc` · `google-cloud-storage`
- GENCODE v41 basic annotation (`/usr/local/gencode.v41.basic.annotation.gff3.gz`), used for ATAC gene-activity conversion
- Scripts: `multiome_label_transfer.py`, `gcs_utils.py` (see [Scripts](#scripts))

## Versioning

This image follows the WARP-Tools tag convention `us.gcr.io/broad-gotc-prod/scvi-scanvi:<image-version>-<scvi-tools-version>-<unix-timestamp>` (see [BUILDING.md → Versioning strategy](../../BUILDING.md#versioning-strategy)). Published tags are recorded in `docker_versions.tsv`. The WARP scANVI pipeline pins the image by immutable `@sha256` digest (shown in [Quick reference](#quick-reference)). Inspect an image with:

```bash
docker inspect us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest>
```

## Scripts

| Script | Language | Purpose |
| --- | --- | --- |
| `multiome_label_transfer.py` | Python | Preprocessing, SCVI/SCANVI training, and label transfer. Exposes importable functions (below) and a CLI `main()`. |
| `gcs_utils.py` | Python | Google Cloud Storage localize/delocalize helpers used by the `--localize` flag. |

Key importable functions in `multiome_label_transfer.py`:

- `run_multi_model(gex, atac_activity, ref)` — train SCVI+SCANVI on GEX + ATAC activity + reference.
- `run_gex_only_model(gex, ref)` — train SCVI+SCANVI on GEX + reference only (GEX-only mode).
- `_concat_and_train_models(adatas, keys)` — shared training core used by both of the above (avoids drift).
- `transfer_labels(data, lvae)` — predict labels (`C_scANVI`) and compute the SCANVI UMAP.
- `finalize_output(scanvi_output)` — add metadata, rename `final_annotation` → `celltype`, copy counts into `.raw`.

## Inputs

The CLI consumes three AnnData `.h5ad` files. (In the WARP pipeline, ATAC is optional — see [Task descriptions](#task-descriptions).)

| Input | Type / format | Required | Description |
| --- | --- | --- | --- |
| `-g` / `--gex-file` | AnnData `.h5ad` | yes | Gene expression. Expected `obs['star_IsCell']` (STARsolo cell calls). |
| `-a` / `--atac-file` | AnnData `.h5ad` (cell-by-bin) | yes (CLI) | ATAC cell-by-bin matrix. Expected `obs['gex_barcodes']` to align with GEX. |
| `-r` / `--ref-file` | AnnData `.h5ad` | yes | Annotated reference. **Must** carry `obs['final_annotation']` (cell-type labels) **and** `obs['batch']` (training is batch-aware). |
| `-l` / `--localize` | flag | no | Localize inputs from / delocalize outputs to GCS via `gcs_utils`. |

## Outputs

| Output | Format | Description |
| --- | --- | --- |
| `SCANVI_predictions.h5ad` | AnnData | Combined object with SCANVI predictions (`celltype`), UMAP, and metadata; counts in `.raw`. |
| `gex_annotated_matrix.h5ad` | AnnData | GEX cells with transferred `final_annotation`. |
| `atac_annotated_matrix.h5ad` | AnnData | ATAC gene-activity cells with transferred `final_annotation` (multiome mode only). |
| `vae_test_model_/`, `scvi_scanvi_test_model_/` | scvi-tools model dirs | Saved SCVI and SCANVI models. |
| `adata_scanvi_labels.h5ad` | AnnData | Intermediate labeled object written by `transfer_labels`. |

## Usage

Run the multiome workflow standalone (GPU required for reasonable runtime):

```bash
docker run --rm -it --gpus all \
    -v /files:/files \
    us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> \
    python multiome_label_transfer.py \
      --gex-file /files/gex.h5ad \
      --atac-file /files/atac.cellbybin.h5ad \
      --ref-file /files/pbmc.h5ad
```

With `podman` and CDI GPU support, replace `--gpus all` with `--device nvidia.com/gpu=all`.

> The CLI `main()` always runs the **multiome** path. **GEX-only** mode is exercised through the WARP pipeline, which imports `run_gex_only_model` directly (it does not call `main()`).

## Task descriptions

This image is consumed by the WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi/scANVI.wdl):

- **`PreprocessFilter`** (CPU) — runs inline Python (loads/filters GEX, builds the ATAC gene-activity matrix with `snapatac2` in multiome mode); does not import from the script.
- **`MultiomeLabelTransfer`** (GPU) — imports `run_multi_model` / `run_gex_only_model`, `transfer_labels`, and `finalize_output` from `/usr/local`, and **never calls `main()`**. It calls `run_multi_model` when ATAC is present and `run_gex_only_model` when it is not.

The imported function API is the production contract — keep signatures and behavior stable (see [AGENTS.md → Keep scripts importable](../../AGENTS.md#keep-scripts-importable)).

## Requirements / dependencies

- **GPU strongly recommended** (e.g. NVIDIA T4 or better); pass `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman/CDI).
- **Network at runtime**: the ATAC gene-activity step downloads `snapatac2`'s hg38 annotation (`snap.genome.hg38`).
- Reference sizing: the PBMC reference is ~2 GB; a full multiome run is comfortable with ~30 GB RAM and ~200 GB disk.

## Troubleshooting

To poke around the image, start a shell: `docker run -it --rm us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> bash`. See [BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
