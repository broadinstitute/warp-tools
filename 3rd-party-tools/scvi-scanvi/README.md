# scvi-scanvi

## Quick reference

```bash
docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:775d02a46203eaecdb12cd99635f1e1a2de20efb37419cb675951ecd56cfb64e
```

- **What is this image:** a GPU-enabled Python image (Debian `python:3.12-slim` base) for single-cell cell-type **label transfer** with deep generative models, bundling the `multiome_label_transfer.py` workflow.
- **What is scvi-tools:** [scvi-tools](https://docs.scvi-tools.org/) provides SCVI (unsupervised variational inference) and SCANVI (semi-supervised annotation) models for single-cell omics.

## Synopsis

This image trains SCVI/SCANVI models to transfer cell-type labels from an annotated reference onto unlabeled query cells. It supports two modes: **multiome** (gene expression + ATAC gene-activity + reference) and **GEX-only** (gene expression + reference, no ATAC). It is the execution environment for the Broad WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi); the workflow follows the [SnapATAC2 annotation tutorial](https://kzhang.org/SnapATAC2/tutorials/annotation.html) and uses the annotated PBMC reference from the scArches tutorial.

## Image contents

- Base image: `python:3.12-slim-trixie` (`--platform=linux/amd64`)
- `scvi-tools` 1.5.1 · `snapatac2` 2.10.0 · `scanpy` · `anndata` · `numpy` · `scikit-misc` · `google-cloud-storage`
- ATAC gene-activity conversion uses snapatac2's hg38 annotation (`snap.genome.hg38`), fetched at runtime. (A GENCODE v41 GFF3 is baked in at `/usr/local/gencode.v41.basic.annotation.gff3.gz` but is legacy — the current code path does not use it.)
- Scripts: `multiome_label_transfer.py`, `label_transfer_from_preprocessed.py`, `gcs_utils.py` (see [Scripts](#scripts))

## Versioning

This image follows the WARP-Tools tag convention `us.gcr.io/broad-gotc-prod/scvi-scanvi:<image-version>-<scvi-tools-version>-<unix-timestamp>` (see [BUILDING.md → Versioning strategy](../../BUILDING.md#versioning-strategy)). Published references are recorded in [docker_versions.tsv](docker_versions.tsv), with the last line being the version currently used by WARP. The current production image is built by CI (branch-tagged) and the scANVI pipeline pins it by immutable `@sha256` digest (shown in [Quick reference](#quick-reference)) — hence the digest rather than a semantic tag for the latest entry. Inspect an image with:

```bash
docker inspect us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest>
```

## Scripts

| Script | Language | Purpose |
| --- | --- | --- |
| `scanvi_sample.py` | Python | **WARP entry point.** Labels one sample end to end: `preprocess_query.py`, then `label_transfer()`. The scANVI WDL runs it once per sample, each in its own process. |
| `validate_inputs.py` | Python | Checks every sample's input files and ids before any sample runs, reporting every problem at once. |
| `preprocess_query.py` | Python | `preprocess_query(gex_path, atac_path, genome)`: one sample's GEX filtering and, in multiome mode, the ATAC gene-activity conversion (formerly the WDL's inline PreprocessFilter task). |
| `label_transfer_from_preprocessed.py` | Python | `label_transfer()`: obtain the model (train, load, or map the sample onto a reference-only model), transfer labels, write outputs. Also a CLI over preprocessed files. `--seed N` makes training reproducible. |
| `scp_output.py` | Python | `make_scp_ready()`: puts an annotated matrix in the form Single Cell Portal reads (see [Outputs](#outputs)). |
| `multiome_label_transfer.py` | Python | SCVI/SCANVI training and label transfer functions (below), plus a standalone CLI `main()`. |
| `make_reference_only_model.py` | Python | Trains a SCANVI model on a reference alone, the way warp-scanvi-training does. Builds the WARP query-mapping test fixture. |
| `test_seed_reproducibility.py` | Python | CI check: two `--seed 42` training runs on synthetic data must give identical labels and embedding. |
| `test_scp_and_query_mapping.py` | Python | CI check: SCP-ready outputs, run provenance and query mapping onto a reference-only model, on synthetic data; and that `validate_inputs.py` rejects bad inputs. |
| `gcs_utils.py` | Python | Google Cloud Storage localize/delocalize helpers used by the `--localize` flag. |

Key importable functions in `multiome_label_transfer.py`:

- `run_multi_model(gex, atac_activity, ref)` — train SCVI+SCANVI on GEX + ATAC activity + reference.
- `run_gex_only_model(gex, ref)` — train SCVI+SCANVI on GEX + reference only (GEX-only mode).
- `_concat_and_train_models(adatas, keys)` — shared training core used by both of the above (avoids drift).
- `transfer_labels(data, lvae)` — predict labels (`C_scANVI`) and compute the SCANVI UMAP.
- `finalize_output(scanvi_output)` — rename `final_annotation` → `celltype`, copy counts into `.raw`. Adds no sample metadata.

## Inputs

The CLI consumes three AnnData `.h5ad` files. (In the WARP pipeline, ATAC is optional — see [Task descriptions](#task-descriptions).)

| Input | Type / format | Required | Description |
| --- | --- | --- | --- |
| `-g` / `--gex-file` | AnnData `.h5ad` | yes | Gene expression. Expected `obs['star_IsCell']` (STARsolo cell calls). |
| `-a` / `--atac-file` | AnnData `.h5ad` (cell-by-bin) | yes (CLI) | ATAC cell-by-bin matrix. Expected `obs['gex_barcodes']` to align with GEX. |
| `-r` / `--ref-file` | AnnData `.h5ad` | yes | Annotated reference. **Must** carry `obs['final_annotation']` (cell-type labels) **and** `obs['batch']` (training is batch-aware). |
| `-e` / `--max-epochs` | int | no | Maximum SCVI/SCANVI training epochs (default `500`). Lower it for fast smoke runs or to bound wall-clock on very large datasets. Also available on `run_multi_model`/`run_gex_only_model` as the `max_epochs` keyword. |
| `-l` / `--localize` | flag | no | Localize inputs from / delocalize outputs to GCS via `gcs_utils`. |

## Outputs

`scanvi_sample.py` and `label_transfer_from_preprocessed.py` write, prefixed by the sample's input id and compressed with gzip:

| Output | Format | Description |
| --- | --- | --- |
| `<id>_gex_annotated_matrix.h5ad` | AnnData | The sample's GEX cells with transferred `final_annotation`, in SCP form (below). |
| `<id>_atac_annotated_matrix.h5ad` | AnnData | The sample's ATAC gene-activity cells, same form (multiome mode only). |
| `<id>_SCANVI_predictions.h5ad` | AnnData | The labelled data the model saw (`celltype`, latent, UMAP; counts in `.raw`): the sample only with a supplied model, the sample plus the reference when training. |
| `<id>_scanvi_model.tar.gz` | tar.gz | The minimal SCANVI model (no bundled data), reusable as `--scanvi-model`. |

SCP form (`make_scp_ready`): counts in `.raw` (float32); `.X` normalized to 10,000 per cell then `log1p` (float32); `obsm['X_scANVI']`, the 30-dimensional latent; `obsm['X_umap']`, 2 columns fitted on the sample's own latent; no `counts` layer. Every output records its run provenance in `uns`: `scanvi_input_id`, `scanvi_input_file` and `scanvi_model`. No sample metadata is added.

The standalone `multiome_label_transfer.py` CLI writes unprefixed `SCANVI_predictions.h5ad`, `gex_annotated_matrix.h5ad`, `atac_annotated_matrix.h5ad`, its model directories and an intermediate `adata_scanvi_labels.h5ad`, in the older form.

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

- **`PreprocessReference`** (CPU, training runs only) — runs inline Python to adapt the reference; does not import from the scripts.
- **`MultiomeLabelTransfer`** (GPU/CPU) — runs `/usr/local/validate_inputs.py` once, then `/usr/local/scanvi_sample.py` once per sample, each in its own process and working directory. A supplied model that knows the sample's batch label is loaded as is; one trained on the reference alone gets the sample mapped onto it (`SCANVI.load_query_data` plus `--query-max-epochs` of training on the sample).

The imported function API is the production contract — keep signatures and behavior stable (see [AGENTS.md → Keep scripts importable](../../AGENTS.md#keep-scripts-importable)).

## Requirements / dependencies

- **GPU strongly recommended** (e.g. NVIDIA T4 or better); pass `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman/CDI).
- **Network at runtime**: the ATAC gene-activity step downloads `snapatac2`'s hg38 annotation (`snap.genome.hg38`).
- Reference sizing: the PBMC reference is ~2 GB; a full multiome run is comfortable with ~30 GB RAM and ~200 GB disk.

## Troubleshooting

To poke around the image, start a shell: `docker run -it --rm us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> bash`. See [BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
