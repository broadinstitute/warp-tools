# scvi-scanvi

## Quick reference

```bash
docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:5764b99116daedef6c73a83483028ac78e318b1ec33de8824a9861b4d3201b52
```

- **What is this image:** a GPU-enabled Python image (Debian `python:3.12-slim` base) for single-cell cell-type **label transfer** with deep generative models, bundling the `multiome_label_transfer.py` workflow.
- **What is scvi-tools:** [scvi-tools](https://docs.scvi-tools.org/) provides SCVI (unsupervised variational inference) and SCANVI (semi-supervised annotation) models for single-cell omics.

## Synopsis

This image trains SCVI/SCANVI models to transfer cell-type labels from an annotated reference onto unlabeled query cells. It supports two modes: **multiome** (gene expression + ATAC gene-activity + reference) and **GEX-only** (gene expression + reference, no ATAC). It is the execution environment for the Broad WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi); the workflow follows the [SnapATAC2 annotation tutorial](https://kzhang.org/SnapATAC2/tutorials/annotation.html) and uses the annotated PBMC reference from the scArches tutorial.

## Image contents

- Base image: `python:3.12-slim-trixie` (`--platform=linux/amd64`)
- `scvi-tools` 1.5.1 and `snapatac2` 2.10.0, pinned; `scanpy`, `anndata`, `numpy`, `torch`, `scipy` and `scikit-misc` at the newest versions compatible with those when the image is built.
- Each build records its exact library versions: `/usr/local/versions.txt` holds the build's `pip freeze`, the CI build log prints it, and every output h5ad records the key versions in `uns['scanvi_versions']`.
- ATAC gene-activity conversion uses snapatac2's genome annotations (`snap.genome.hg38`, `mm10`, `mm39`), downloaded at runtime.
- Scripts: see [Scripts](#scripts).

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
| `reference_id.py` | Python | `reference_id()`: names a reference by its `uns['title']` (or file name) plus the first 8 hex digits of the file's sha256; models carry it in `reference.txt`. |
| `scp_output.py` | Python | `make_scp_ready()`: puts an annotated matrix in the form Single Cell Portal reads, with `.X` as log2(CPM + 1) (see [Outputs](#outputs)). |
| `multiome_label_transfer.py` | Python | SCVI/SCANVI training and label transfer functions (below), imported by `label_transfer_from_preprocessed.py`. No CLI. |
| `make_reference_only_model.py` | Python | Trains a SCANVI model on a reference alone, the way warp-scanvi-training does. Builds the WARP query-mapping test fixture. |
| `test_seed_reproducibility.py` | Python | CI check: two `--seed 42` training runs on synthetic data must give identical labels and embedding. |
| `test_scp_and_query_mapping.py` | Python | CI check: SCP-ready outputs, run provenance and query mapping onto a reference-only model, on synthetic data; and that `validate_inputs.py` rejects bad inputs. |

Key importable functions in `multiome_label_transfer.py`:

- `run_multi_model(gex, atac_activity, ref)` — train SCVI+SCANVI on GEX + ATAC activity + reference.
- `run_gex_only_model(gex, ref)` — train SCVI+SCANVI on GEX + reference only (GEX-only mode).
- `_concat_and_train_models(adatas, keys)` — shared training core used by both of the above (avoids drift).
- `transfer_labels(data, lvae)` — predict labels (`C_scANVI`), store the latent (`X_scANVI`), and fit neighbors and a UMAP on it. Returns `data`.
- `finalize_output(scanvi_output)` — rename `final_annotation` → `celltype`, copy counts into `.raw`. Adds no sample metadata.

## Inputs

`scanvi_sample.py` labels one sample:

| Argument | Required | Description |
| --- | --- | --- |
| `--input-id` | yes | The sample's id; it prefixes every output file and is recorded in each output. |
| `--gex` | yes | The sample's GEX h5ad (Optimus output). A missing `obs['star_IsCell']` is added as all True. |
| `--atac` | no | The sample's ATAC cell-by-bin h5ad (multiome mode); barcodes are matched through `obs['gex_barcodes']` (or the index). |
| `--genome` | no | Genome for the ATAC gene-activity conversion: `hg38` (default), `mm10` or `mm39`. |
| `--scanvi-model` | one of these two | A saved SCANVI model (`.tar.gz`), used instead of training. |
| `--ref` | one of these two | A reference preprocessed by the WDL's PreprocessReference task (counts in `.X`, labels in `obs['final_annotation']`, batch in `obs['batch']`), to train on. |
| `--max-epochs`, `--query-max-epochs`, `--batch-size`, `--seed`, `--output-max-probability` | no | Training epochs (default 500), query-mapping epochs (default 100), minibatch size (default 128), random seed, and a `max_probability` column. |

## Outputs

`scanvi_sample.py` and `label_transfer_from_preprocessed.py` write, prefixed by the sample's input id and compressed with gzip:

| Output | Format | Description |
| --- | --- | --- |
| `<id>_gex_annotated_matrix.h5ad` | AnnData | The sample's GEX cells with transferred `final_annotation`, in SCP form (below). |
| `<id>_atac_annotated_matrix.h5ad` | AnnData | The sample's ATAC gene-activity cells, same form (multiome mode only). |
| `<id>_SCANVI_predictions.h5ad` | AnnData | The labelled data the model saw (`celltype`, latent, UMAP; counts in `.raw`): the sample only with a supplied model, the sample plus the reference when training. |
| `<id>_scanvi_model.tar.gz` | tar.gz | The minimal SCANVI model (no bundled data) plus `reference.txt`, reusable as `--scanvi-model`. |

SCP form (`make_scp_ready`): counts in `.raw` (float32); `.X` is log2(CPM + 1) (float32; units are bits; Allen's convention, `uns['log1p']['base'] = 2`); `obsm['X_scANVI']`, the 30-dimensional latent; `obsm['X_umap']`, 2 columns fitted on the sample's own latent; no `counts` layer. Every output records its run provenance in `uns`: `scanvi_input_id`, `scanvi_input_file` and `scanvi_model`, the model's name, which is its reference's id from `reference.txt` (or the archive's file name when it has none). No sample metadata is added.

## Usage

Label one sample with a saved model, on CPU:

```bash
docker run --rm -v /files:/files -w /files \
    us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> \
    python3 /usr/local/scanvi_sample.py \
      --input-id sample_a \
      --gex /files/sample_a_gex.h5ad \
      --scanvi-model /files/scanvi_model.tar.gz
```

To train instead, pass `--ref` with a reference prepared by the WDL's PreprocessReference task, and attach a GPU: `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman with CDI). The WARP pipeline is the normal way to run this image.

## Task descriptions

This image is consumed by the WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi/scANVI.wdl):

- **`PreprocessReference`** (CPU, training runs only) — runs inline Python to adapt the reference; does not import from the scripts.
- **`MultiomeLabelTransfer`** (GPU/CPU) — runs `/usr/local/validate_inputs.py` once, then `/usr/local/scanvi_sample.py` once per sample, each in its own process and working directory. A supplied model that knows the sample's batch label is loaded as is; one trained on the reference alone gets the sample mapped onto it (`SCANVI.load_query_data` plus `--query-max-epochs` of training on the sample).

The imported function API is the production contract — keep signatures and behavior stable (see [AGENTS.md → Keep scripts importable](../../AGENTS.md#keep-scripts-importable)).

## Requirements / dependencies

- **GPU strongly recommended** (e.g. NVIDIA T4 or better); pass `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman/CDI).
- **Network at runtime**: the ATAC gene-activity step downloads `snapatac2`'s genome annotation.
- Reference sizing: the PBMC reference is ~2 GB; a full multiome run is comfortable with ~30 GB RAM and ~200 GB disk.

## Troubleshooting

To poke around the image, start a shell: `docker run -it --rm us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> bash`. See [BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
