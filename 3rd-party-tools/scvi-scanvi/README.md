# Label Transfer for Multiome data with SCVI and SCANVI
This image outlines a workflow for integrating single-cell RNA and ATAC data using deep generative models (SCVI and SCANVI).
It demonstrates how to preprocess the data, train models, and transfer cell type labels from a gene expression (GEX) 
reference dataset to an ATAC-seq query dataset for downstream analysis and visualization.

We use a full Multiome dataset with peak calling as our query, and the annotated PBMC reference from the scArches tutorial.

This workflow follows the SnapATAC2 annotation tutorial: SnapATAC2 Tutorial https://kzhang.org/SnapATAC2/tutorials/annotation.html
## Environment Setup
Environment Setup
Note: GPU acceleration is required to run this tutorial efficiently. Ensure that GPU support is enabled when creating your environment.

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