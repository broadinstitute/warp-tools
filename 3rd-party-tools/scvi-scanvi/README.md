# Label Transfer for Multiome data with SCVI and SCANVI
This image outlines a workflow for integrating single-cell RNA and ATAC data using deep generative models (SCVI and SCANVI).
It demonstrates how to preprocess the data, train models, and transfer cell type labels from a gene expression (GEX) 
reference dataset to an ATAC-seq query dataset for downstream analysis and visualization.

We use a full Multiome dataset with peak calling as our query, and the annotated PBMC reference from the scArches tutorial.

This workflow follows the SnapATAC2 annotation tutorial: SnapATAC2 Tutorial https://kzhang.org/SnapATAC2/tutorials/annotation.html

## Quick reference

Copy and paste to pull this image

#### `docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975`

## Versioning

This image uses the following convention for versioning:

#### `us.gcr.io/broad-gotc-prod/scvi-scanvi:<image-version>-<scvi-version>-<unix-timestamp>` 

We keep track of all past versions in [docker_versions](docker_versions.tsv) with the last image listed being the currently used version in WARP.

You can see more information about the image, including the tool versions, by running the following commands:

```bash
$ docker pull us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975
$ docker inspect us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975
```

## Usage
```bash
$ docker run --rm -it \
    -v /files:/files --gpus all \
    us.gcr.io/broad-gotc-prod/scvi-scanvi:1.0.0-1.2-1756234975 \
    python multiome_label_transfer.py --gex-file /files/gex.h5ad \
    --atac-file /files/atac.h5ad --ref-file /files/pbmc.h5ad 
```


## Environment Setup
Environment Setup
Note: GPU acceleration is strongly recommended in order to run this container efficiently, though it is not a strict hardware requirement. 
Ensure that GPU support is enabled when creating your environment. When running the container, you must grant access to 
any GPUs by setting the `--gpus all` flag in `docker run`.

VM Specifications:

* GPU: NVIDIA Tesla T4
* CPU: 8 vCPUs
* Memory: 30 GB RAM
* Disk: 200 GB
* Input Files and Sizes:
  * GEX h5ad file test_gex.h5ad: 675 MB
  * ATAC peaks called per cluster human.cellbybin.h5ad: 591 MB 
  * Annotated reference pbmc.h5ad: 2 GB

## Overview

### Train SCVI and SCANVI models to integrate RNA (query and reference) and ATAC datasets
#### Preprocessing steps:

Filter gene expression (GEX) using `filt_gex` to retain barcodes called as cells by `STARsolo` and store raw counts in .layers.
Reindex ATAC using `reindex_barcodes` to align ATAC data with GEX barcodes to ensure matching cells.
Filter to shared barcodes using `shared_barcodes` to retain only the intersecting cells in GEX and ATAC datasets.
Align gene features using `get_shared_features` to ensure GEX, ATAC (activity matrix), and reference datasets share the same gene space.
#### Training process:

Combine GEX, ATAC, and reference & select highly variable genes using `run_multi_model`.
Train SCVI model using `run_multi_model` to learn an unsupervised latent representation from the combined data.
Train SCANVI model using `run_multi_model` to transfer cell type labels from reference to query using semi-supervised learning.
Save outputs of trained models and annotated AnnData objects.
#### Transfer labels using SCANVI model.
Use transfer_labels to apply the trained SCANVI model to new data and visualize cell types via UMAP.

### Explore previously trained model for Single Cell Portal.
Separate the combined multiome AnnData object into RNA and ATAC subsets (rna_unannotated and atac_unannotated), while preserving:

transferred annotations
shared barcodes
relevant metadata for downstream analysis.
### Perform label transfer and save predictions.
Use:

`subset_multiome_data` to isolate modalities
`check_cell_matches` to validate cell integrity across modalities.
### Output timing summary table
Report time taken for each stage: preprocessing, training, and prediction.