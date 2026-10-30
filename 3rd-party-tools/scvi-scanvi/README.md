# scvi-scanvi

## Quick reference

```bash
docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:938baa80261c9d2122838403d36e000cf914035705f4cbddc6eed34595289386
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

This image follows the WARP-Tools tag convention `us.gcr.io/broad-gotc-prod/scvi-scanvi:<image-version>-<scvi-tools-version>-<unix-timestamp>` (see [BUILDING.md → Versioning strategy](../../BUILDING.md#versioning-strategy)). Published references are recorded in [docker_versions.tsv](docker_versions.tsv), with the last line being the version currently used by WARP. What changed in each image version, and the library versions of the current one, are in [scvi-scanvi.changelog.md](scvi-scanvi.changelog.md). The current production image is built by CI (branch-tagged) and the scANVI pipeline pins it by immutable `@sha256` digest (shown in [Quick reference](#quick-reference)) — hence the digest rather than a semantic tag for the latest entry. Inspect an image with:

```bash
docker inspect us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest>
```

## Scripts

| Script | Language | Purpose |
| --- | --- | --- |
| `scanvi_sample.py` | Python | **WARP entry point.** Labels one sample end to end: `preprocess_query.py`, then `label_transfer()`. The scANVI WDL runs it once per sample, each in its own process. |
| `validate_inputs.py` | Python | Checks every sample's input files and ids, and that no two outputs share a name, before any sample runs, reporting every problem at once. |
| `preprocess_query.py` | Python | `preprocess_query(gex_path, atac_path, genome)`: one sample's GEX filtering and, in multiome mode, the ATAC gene-activity conversion (formerly the WDL's inline PreprocessFilter task). |
| `label_transfer_from_preprocessed.py` | Python | `label_transfer()`: for each label level, obtain the model (train, load, or map the sample onto a reference-only model) and transfer labels; then write outputs. Also a CLI over preprocessed files. `--seed N` makes training reproducible. |
| `model_archive.py` | Python | The model archive (below): `write_archive()`, `read_archive()`, and `reference_id()`, which names a reference by its `uns['title']` (or file name) plus the first 8 hex digits of the file's sha256. warp-scanvi-training keeps a copy. |
| `scp_output.py` | Python | `make_scp_ready()`: puts an annotated matrix in the form Single Cell Portal reads, with `.X` as log2(CPM + 1) (see [Outputs](#outputs)). |
| `multiome_label_transfer.py` | Python | SCVI/SCANVI training on a reference plus one sample (below), imported by `label_transfer_from_preprocessed.py`. No CLI. |
| `make_reference_only_model.py` | Python | Trains scVI and a SCANVI per label level on a reference alone, the way warp-scanvi-training does, into a model archive. Builds the WARP query-mapping test fixture. |
| `test_seed_reproducibility.py` | Python | CI check: two `--seed 42` training runs on synthetic data must give identical labels and embedding. |
| `test_scp_and_query_mapping.py` | Python | CI check: SCP-ready outputs, run provenance and query mapping onto a two-level reference-only model, on synthetic data; and that `validate_inputs.py` rejects bad inputs. |

Key importable functions in `multiome_label_transfer.py`:

- `train_scvi(gex, atac_activity, ref)` — concatenate the sample (GEX, plus ATAC activity in multiome mode) with the reference and train SCVI. Returns `(data, vae)`.
- `train_scanvi(data, vae, ref, level)` — train a SCANVI from `vae` with the reference's `obs[level]` as labels; the sample's cells are unlabelled.

### Model archive

One `.tar.gz` holds every model trained from one reference (`model_archive.py`):

```text
scanvi_model/
  reference.txt             the reference id: the model's name in every output
  levels.txt                the label levels, one per line, in order; the first is SCP's default
  scvi/model.pt             the scVI model the SCANVI models were trained from (optional)
  scanvi/<level>/model.pt   one SCANVI model per label level, e.g. class, subclass
```

Models are saved without data. Level names are limited to letters, digits, `_`, `.` and `-`, because each becomes an `obs` column and part of `obsm` keys.

## Inputs

`scanvi_sample.py` labels one sample:

| Argument | Required | Description |
| --- | --- | --- |
| `--input-id` | yes | The sample's id, recorded in each output. |
| `--gex` | yes | The sample's GEX h5ad (Optimus output). A missing `obs['star_IsCell']` is added as all True. |
| `--gex-out` | yes | Path of the annotated GEX h5ad to write (the WDL uses `<gex file name>_annotated.h5ad`). |
| `--atac`, `--atac-out` | no | The sample's ATAC cell-by-bin h5ad (multiome mode; barcodes are matched through `obs['gex_barcodes']` or the index), and the annotated ATAC h5ad to write. |
| `--genome` | no | Genome for the ATAC gene-activity conversion: `hg38` (default), `mm10` or `mm39`. |
| `--scanvi-model` | one of these two | A model archive (`.tar.gz`), used instead of training. A model whose genes are Ensembl IDs is matched through the sample's `ensembl_ids` (Optimus) column; the log reports how many of the model's genes the sample has. |
| `--ref` | one of these two | A reference preprocessed by the WDL's PreprocessReference task (counts in `.X`, one string column per label level listed in `uns['scanvi_levels']`, batch in `obs['batch']`), to train on. |
| `--model-out` | no | When training, also write the trained models as a model archive. |
| `--max-epochs`, `--query-max-epochs`, `--batch-size`, `--seed`, `--output-max-probability` | no | Training epochs (default 500), query-mapping epochs (default 100), minibatch size (default 128), random seed, and a `<level>_max_probability` column per level. |

## Outputs

`scanvi_sample.py` and `label_transfer_from_preprocessed.py` write to the paths given, compressed with gzip:

| Output | Format | Description |
| --- | --- | --- |
| `--gex-out` | AnnData | The sample's GEX cells with their labels at each level, in SCP form (below). |
| `--atac-out` | AnnData | The sample's ATAC gene-activity cells, same form (multiome mode only). |
| `--model-out` | tar.gz | Training only: the model archive (scVI plus a SCANVI per level, no bundled data), reusable as `--scanvi-model`. |

SCP form (`make_scp_ready`): counts in `.raw` (float32); `.X` is log2(CPM + 1) (float32; units are bits; Allen's convention, `uns['log1p']['base'] = 2`); per level, the labels in `obs['<level>']`, the 30-dimensional scVI latent in `obsm['X_scVI_<level>']` and a 2-column UMAP fitted on it, for the sample's own cells, in `obsm['X_umap_<level>']`; `obsm['X_umap']` is a copy of the first level's; no `counts` layer. Every output records its run provenance in `uns`: `scanvi_input_id`, `scanvi_input_file`, `scanvi_model` (the model's name, its reference's id from `reference.txt`), `scanvi_levels` and `scanvi_versions`. No sample metadata is added.

## Usage

Label one sample with a saved model, on CPU:

```bash
docker run --rm -v /files:/files -w /files \
    us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> \
    python3 /usr/local/scanvi_sample.py \
      --input-id sample_a \
      --gex /files/sample_a_gex.h5ad \
      --gex-out /files/sample_a_gex_annotated.h5ad \
      --scanvi-model /files/scanvi_model.tar.gz
```

To train instead, pass `--ref` with a reference prepared by the WDL's PreprocessReference task, and attach a GPU: `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman with CDI). The WARP pipeline is the normal way to run this image.

## Task descriptions

This image is consumed by the WARP [scANVI pipeline](https://github.com/broadinstitute/warp/tree/develop/pipelines/wdl/scanvi/scANVI.wdl):

- **`PreprocessReference`** (CPU, training runs only) — runs inline Python to adapt the reference; imports `check_levels()` and `reference_id()` from `/usr/local/model_archive.py`.
- **`MultiomeLabelTransfer`** (GPU/CPU) — runs `/usr/local/validate_inputs.py` once, then `/usr/local/scanvi_sample.py` once per sample, each in its own process and working directory. For each level, a supplied model that knows the sample's batch label is loaded as is; one trained on the reference alone gets the sample mapped onto it (`SCANVI.load_query_data` plus `--query-max-epochs` of training on the sample).

The imported function API is the production contract — keep signatures and behavior stable (see [AGENTS.md → Keep scripts importable](../../AGENTS.md#keep-scripts-importable)).

## Requirements / dependencies

- **GPU strongly recommended** (e.g. NVIDIA T4 or better); pass `--gpus all` (docker) or `--device nvidia.com/gpu=all` (podman/CDI).
- **Network at runtime**: the ATAC gene-activity step downloads `snapatac2`'s genome annotation.
- Reference sizing: the PBMC reference is ~2 GB; a full multiome run is comfortable with ~30 GB RAM and ~200 GB disk.

## Troubleshooting

To poke around the image, start a shell: `docker run -it --rm us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:<digest> bash`. See [BUILDING.md → Troubleshooting and running standalone](../../BUILDING.md#troubleshooting-and-running-standalone).
