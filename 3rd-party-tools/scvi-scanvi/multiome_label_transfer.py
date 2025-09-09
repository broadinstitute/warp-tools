import warnings
warnings.filterwarnings("ignore")

import argparse
import snapatac2 as snap
import anndata as ad
import numpy as np
import scanpy as sc
import scvi
import time
from gcs_utils import *

def read_h5ad(gex_path, cellbybin_path, reference_path):
    """
        Load gene expression, chromatin accessibility, and reference datasets.

        Parameters:
        - gex_path (str): Path to the gene expression dataset (.h5ad).
        - cellbybin_path (str): Path to the ATAC-seq cell-by-bin dataset (Snap format).
        - reference_path (str): Path to the reference dataset (.h5ad).

        Returns:
        - gex (AnnData): Gene expression data.
        - atac (SnapATAC object): Chromatin accessibility data.
        - ref (AnnData): Reference dataset.
    """

    gex = sc.read_h5ad(gex_path)
    atac = snap.read((cellbybin_path), backed=None)
    ref = sc.read_h5ad(reference_path)
    print("GEX adata: ", gex.shape)
    print("ATAC adata: ", atac.shape)
    print("Refence adata: ", ref.shape)
    return gex, atac, ref


def filt_gex(gex_adata):
    """
        Filter a gene expression AnnData object to retain high-quality cells and genes.
        Filters Multiome output GEX h5ad so that it only uses barcodes that have been called as cell by STARsolo.
        Adds the counts into adatalayer.

        Steps:
        - Keep only cells marked as true STAR calls (gex.obs['star_IsCell'] == True).
        - Filter out genes and cells with fewer than 3 total counts.
        - Add a 'batch' column with a constant value (1).
        - Copy the expression matrix to a new layer called 'counts'.
        - Save the filtered object to 'filtered_gex.h5ad'.

        Parameters:
        - gex_adata (AnnData): Raw gene expression AnnData object.

        Returns:
        - gex (AnnData): Filtered and preprocessed AnnData object.
    """

    gex = gex_adata
    print("gex shape", gex.shape)

    print("Filtering to STAR cell calls")
    gex = gex[gex.obs['star_IsCell'] == True]

    print("# cells, # genes before filtering:", gex.shape)

    # Filter genes and cells with fewer than 3 counts
    sc.pp.filter_genes(gex, min_counts=3)
    sc.pp.filter_cells(gex, min_counts=3)

    print("# cells, # genes after filtering:", gex.shape)

    # Add batch info and raw count layer
    gex.obs["batch"] = 1
    gex.layers["counts"] = gex.X.copy()

    # Save filtered object
    gex.write("filtered_gex.h5ad")

    return (gex)

def reindex_barcodes(cellbybin_adata):
    """
        Reindex the ATAC AnnData object using gene expression (GEX) barcodes.

        This is typically done to align ATAC cells with corresponding GEX cells
        in Multiome datasets where matched modalities are processed separately.

        Parameters:
        - cellbybin_adata (AnnData): Chromatin accessibility AnnData object (e.g., cell-by-bin matrix).

        Returns:
        - atac (AnnData): The same AnnData object with its obs index set to 'gex_barcodes'.
    """
    atac = cellbybin_adata
    atac.obs = atac.obs.set_index("gex_barcodes")
    return atac

def shared_barcodes(rna, atac):
    """
        Filter ATAC and GEX AnnData objects to retain only cells (barcodes) shared across both modalities.

        This ensures that downstream multimodal analysis is performed on matched cells.

        Parameters:
        - rna (AnnData): Gene expression data.
        - atac (AnnData): Chromatin accessibility (e.g., cell-by-bin) data.

        Returns:
        - gex_shared (AnnData): Filtered GEX AnnData object with only shared barcodes.
        - atac_shared (AnnData): Filtered ATAC AnnData object with only shared barcodes.

        Side Effects:
        - Writes the filtered objects to "gex_filtered.h5ad" and "atac_filtered.h5ad".
    """
    # Find the shared barcodes
    shared_barcodes = atac.obs.index.intersection(rna.obs.index)

    # Subset both AnnData objects
    atac_shared = atac[atac.obs.index.isin(shared_barcodes)].copy()
    gex_shared = rna[rna.obs.index.isin(shared_barcodes)].copy()

    # Save the filtered AnnData objects
    atac_shared.write_h5ad("atac_filtered.h5ad")
    gex_shared.write_h5ad("gex_filtered.h5ad")

    print(f"Filtered ATAC and GEX files saved with {len(shared_barcodes)} shared barcodes.")
    return gex_shared, atac_shared


def get_shared_features(gex, atac, ref):
    """
        Align gene features across GEX, ATAC, and reference AnnData objects.

        This function ensures that all datasets share a consistent gene space for
        downstream multi-modal integration or analysis.

        Assumes that the ATAC activity matrix includes gene features (not peaks/bins).

        Parameters:
        - gex (AnnData): Gene expression data.
        - atac (AnnData): ATAC data (e.g., gene activity matrix).
        - ref (AnnData): Reference dataset.

        Returns:
        - gex (AnnData): Subsetted GEX data with shared gene features.
        - atac (AnnData): Subsetted ATAC data with shared gene features.
        - ref (AnnData): Subsetted reference data with shared gene features.
    """
    print("Aligning gene and region features...")

    # Find gene features in RNA datasets
    gene_features = ref.var_names.intersection(gex.var_names)

    # Find overlapping gene features in ATAC dataset
    atac_features = atac.var_names.intersection(gene_features)

    print("atac_features length", len(atac_features))
    print("gene_features length", len(gene_features))

    # Ensure all datasets share the same gene space
    gex = gex[:, gene_features]
    ref = ref[:, gene_features]
    atac = atac[:, atac_features]  # ATAC activity matrix does contain gene features

    return gex, atac, ref


def run_multi_model(gex, atac_activity_matrix, ref):
    """
        Train a SCVI/SCANVI model on concatenated GEX, ATAC, and reference datasets.

        This function:
        - Merges GEX, ATAC activity, and reference datasets into a unified AnnData object.
        - Identifies highly variable genes (HVGs) for feature selection.
        - Trains a SCVI model for unsupervised representation learning.
        - Transfers cell type annotations using SCANVI (semi-supervised learning).
        - Saves both trained models and annotated data.

        Parameters:
        - gex (AnnData): Gene expression data (unannotated).
        - atac_activity_matrix (AnnData): ATAC data in gene activity format (unannotated).
        - ref (AnnData): Annotated reference GEX dataset with `final_annotation` in `.obs`.

        Returns:
        - data (AnnData): Combined and annotated data object.
        - vae (SCVI model): Trained SCVI model.
        - lvae (SCANVI model): Trained SCANVI model with transferred labels.
    """
    # Concatenate the datasets across modalities
    data = ad.concat(
        [gex, atac_activity_matrix, ref],
        join='inner',
        label='modality',
        keys=["rna_unannotated", "atac_unannotated", "rna_annotated"],
        index_unique='_',
    )
    # Filter genes expressed in at least 5 cells
    sc.pp.filter_genes(data, min_cells=5)

    # Select highly variable genes using the Seurat v3 method
    sc.pp.highly_variable_genes(
        data,
        n_top_genes=5000,
        flavor="seurat_v3",
        batch_key="batch",
        subset=True
    )

    # Setup SCVI model
    scvi.model.SCVI.setup_anndata(data, batch_key="batch")
    vae = scvi.model.SCVI(
        data,
        n_layers=2,
        n_latent=30,
        gene_likelihood="nb",
        dispersion="gene-batch",
    )

    # Train SCVI
    vae.train(max_epochs=500, early_stopping=True)
    vae.save("vae_test_model_", save_anndata=True)

    # Plot training history
    ax = vae.history['elbo_train'][1:].plot()
    vae.history['elbo_validation'].plot(ax=ax)

    # Prepare labels for SCANVI
    data.obs["celltype_scanvi"] = 'Unknown'
    ref_idx = data.obs['modality'] == "rna_annotated"
    data.obs["celltype_scanvi"][ref_idx] = data.obs['final_annotation'][ref_idx]

    # Initialize and train SCANVI
    lvae = scvi.model.SCANVI.from_scvi_model(
        vae,
        adata=data,
        labels_key="celltype_scanvi",
        unlabeled_category="Unknown",
    )
    lvae.train(max_epochs=500, n_samples_per_label=100)
    lvae.save("scvi_scanvi_test_model_", save_anndata=True)

    return data, vae, lvae


def transfer_labels(data, lvae):
    """
        Use a trained SCANVI model to transfer cell type labels and generate UMAP embedding.

        Assumes that `lvae` (a trained SCANVI model) is already defined in the environment.

        Steps:
        - Predict cell types using SCANVI.
        - Extract SCANVI latent representation.
        - Compute neighbors and UMAP embedding using the SCANVI latent space.
        - Save labeled and embedded AnnData to disk.
        - Return the UMAP plot and the updated AnnData object.

        Parameters:
        - data (AnnData): Combined dataset passed to SCANVI model (from run_multi_model).

        Returns:
        - plot: The UMAP plot object.
        - data (AnnData): Updated AnnData with predicted labels and SCANVI embedding.
    """
    # Predict labels using SCANVI
    data.obs["C_scANVI"] = lvae.predict(data)

    # Get latent representation from SCANVI
    data.obsm["X_scANVI"] = lvae.get_latent_representation(data)

    # Compute neighborhood graph and UMAP based on SCANVI latent space
    sc.pp.neighbors(data, use_rep="X_scANVI")
    sc.tl.umap(data)

    # Plot UMAP colored by predicted labels and modality
    plot = sc.pl.umap(data, color=['C_scANVI', "modality"], wspace=0.45)

    # Save the updated AnnData object
    data.write("adata_scanvi_labels.h5ad")

    return plot, data


def subset_multiome_data(scanvi_output):
    """
        Extracts and prepares RNA and ATAC unannotated subsets from a SCANVI-labeled AnnData object.

        This function:
        - Copies SCANVI-predicted labels to a new field.
        - Separates RNA and ATAC unannotated modalities.
        - Strips modality suffixes from barcodes.
        - Filters to shared barcodes across RNA and ATAC.
        - Renames labels, adds metadata fields for downstream compatibility.
        - Stores raw count layers.
        - Saves processed AnnData objects to disk.

        Parameters:
        - scanvi_output (AnnData): SCANVI output containing multiple modalities and predicted labels.

        Returns:
        - rna_unannotated (AnnData): Processed RNA unannotated dataset with metadata and raw counts.
        - atac_unannotated (AnnData): Processed ATAC unannotated dataset with metadata and raw counts.
    """

    # Copy predicted SCANVI labels to a standard field
    adata = scanvi_output
    adata.obs['final_annotation'] = adata.obs['C_scANVI'].to_numpy()

    # Subset unannotated RNA and ATAC data
    rna_unannotated = adata[adata.obs["modality"] == "rna_unannotated"].copy()
    atac_unannotated = adata[adata.obs["modality"] == "atac_unannotated"].copy()

    # Remove modality suffixes from barcodes
    rna_unannotated.obs.index = rna_unannotated.obs.index.str.replace("_rna_unannotated", "", regex=False)
    atac_unannotated.obs.index = atac_unannotated.obs.index.str.replace("_atac_unannotated", "", regex=False)

    # Keep only shared barcodes between RNA and ATAC
    common_barcodes = rna_unannotated.obs.index.intersection(atac_unannotated.obs.index)
    rna_unannotated = rna_unannotated[common_barcodes]
    atac_unannotated = atac_unannotated[common_barcodes]

    # Save interim files
    rna_unannotated.write_h5ad("rna_unannotated.h5ad")
    atac_unannotated.write_h5ad("atac_unannotated.h5ad")

    # Rename prediction field and add empty ontology placeholders
    rna_unannotated.obs.rename(columns={"final_annotation": "celltype"}, inplace=True)
    atac_unannotated.obs.rename(columns={"final_annotation": "celltype"}, inplace=True)

    # Add gex required metadata fields
    gex_metadata_defaults = {
        "biosample_id": "sample_001",
        "donor_id": "donor_001",
        "species": "NCBITaxon_9606",
        "species__ontology_label": "Homo sapiens",
        "disease": "PATO_0000461",
        "disease__ontology_label": "normal",
        "organ": "UBERON_0000955",
        "organ__ontology_label": "brain",
        "library_preparation_protocol": "EFO_0030059",
        "library_preparation_protocol__ontology_label": "Simultaneous profiling of gene expression and open chromatin from the same cell.",
        "sex": "unknown"
    }

    # Add atac required metadata fields
    atac_metadata_defaults = {
        "biosample_id": "sample_001",
        "donor_id": "donor_001",
        "species": "NCBITaxon_9606",
        "species__ontology_label": "Homo sapiens",
        "disease": "PATO_0000461",
        "disease__ontology_label": "normal",
        "organ": "UBERON_0000955",
        "organ__ontology_label": "brain",
        "library_preparation_protocol": "EFO_0030059",
        "library_preparation_protocol__ontology_label": "Simultaneous profiling of gene expression and open chromatin from the same cell.",
        "sex": "unknown"
    }

    print("Adding metadata to rna h5ad")
    for key, value in gex_metadata_defaults.items():
        rna_unannotated.obs[key] = value

    print("Adding metadata to atac h5ad")
    for key, value in atac_metadata_defaults.items():
        atac_unannotated.obs[key] = value

    # Move raw gene counts and gene activity counts into raw
    rna_unannotated.raw = ad.AnnData(rna_unannotated.X, var=rna_unannotated.var, obs=rna_unannotated.obs)
    atac_unannotated.raw = ad.AnnData(atac_unannotated.X, var=atac_unannotated.var, obs=atac_unannotated.obs)

    # Save the separate AnnData objects
    rna_unannotated.write_h5ad("rna_unannotated_3.h5ad")
    atac_unannotated.write_h5ad("atac_unannotated_3.h5ad")

    return rna_unannotated, atac_unannotated


def check_cell_matches(gex, atac):
    """
        Compare cell type annotations between matched GEX and ATAC AnnData objects
        from the same Multiome library.

        This function assumes that both `gex` and `atac` objects have:
        - Matching cell barcodes (after alignment)
        - A 'celltype' column in .obs containing predicted cell types

        Parameters:
        - gex (AnnData): RNA (gene expression) AnnData object
        - atac (AnnData): ATAC AnnData object

        Returns:
        - mismatched_barcodes (Index): Pandas Index of barcodes with mismatched annotations
    """
    # Find shared barcodes between RNA and ATAC
    # Ensure barcodes match between the two datasets
    common_barcodes = gex.obs.index.intersection(atac.obs.index)

    # Check for mismatches in cell type annotations
    mismatched_barcodes = common_barcodes[
        atac.obs["celltype"].astype(str) != gex.obs["celltype"].astype(str)]

    if not mismatched_barcodes.empty:
        print("Mismatched cell types found for the following barcodes:")
        print(mismatched_barcodes.tolist())
    else:
        print("All cell types match between RNA and ATAC datasets.")

    return mismatched_barcodes


def main(gex_file, atac_file, ref_file):
    # This block performs initial preprocessing steps for integrating Multiome RNA and ATAC data
    # with a reference scRNA-seq dataset using scvi-tools.
    #
    # - Reads in gene expression (GEX), ATAC (cell-by-bin), and reference h5ad files.
    # - Filters GEX cells based on STARsolo cell calls and minimum counts.
    # - Reindexes ATAC cells to match GEX barcodes.
    # - Subsets the data to only the shared barcodes between GEX and ATAC.
    # - Adds metadata such as 'batch' and 'modality' required for downstream SCVI/SCANVI integration.
    # - Converts the ATAC cell-by-bin matrix into a gene activity matrix using the hg38 genome.
    # - Ensures all datasets include a 'final_annotation' field for SCANVI training.

    # Begin loading and filtering...
    timing_summary = {}
    start = time.time()

    # Load the GEX, ATAC, and reference datasets into memory
    gex, atac, ref = read_h5ad(gex_file, atac_file, ref_file)

    # Filter the GEX data to retain only high-quality cells and genes
    gex_filt = filt_gex(gex)
    print("GEX_filtered: ", gex_filt.shape)

    # Reindex the ATAC data so that its barcodes match the GEX barcodes (important for alignment)
    atac = reindex_barcodes(atac)
    print("ATAC n_cells, n_features", atac.shape)
    print("RNA n_cells, n_features", gex_filt.shape)

    # Subset both GEX and ATAC to only the barcodes they have in common (shared cells)
    gex_shared, atac_shared, = shared_barcodes(gex_filt, atac)
    print("GEX shared barcodes: ", gex_shared.shape)
    print("ATAC shared barcodes: ", atac_shared.shape)

    # Assign batch labels for downstream integration with scVI/scanVI
    atac_shared.obs['batch'] = "pd-multiome_sci_atac"
    gex_shared.obs['batch'] = "pd-multiome_sci_gex"

    # Add placeholder annotations to the query datasets (used later by scanVI)
    gex_shared.obs['final_annotation'] = "Unknown"

    query = snap.pp.make_gene_matrix(
        atac_shared,
        gene_anno="gencode.v41.basic.annotation.gff3.gz"
    )

    print("query shape: ", query.shape)
    query.obs['final_annotation'] = "Unknown"

    # Tag each dataset with its modality (required for scanVI integration)
    gex_shared.obs["modality"] = "rna_unannotated"
    query.obs["modality"] = "atac_unannotated"
    ref.obs["modality"] = "rna_annotated"

    timing_summary['Preprocessing'] = time.time() - start

    # Convert ATAC cell-by-bin matrix into gene activity matrix using reference genome annotation
    query = snap.pp.make_gene_matrix(atac_shared, gene_anno=snap.genome.hg38)
    print("query shape: ", query.shape)
    query.obs['final_annotation'] = "Unknown"

    # Tag each dataset with its modality (required for scanVI integration)
    gex_shared.obs["modality"] = "rna_unannotated"
    query.obs["modality"] = "atac_unannotated"
    ref.obs["modality"] = "rna_annotated"

    timing_summary['Preprocessing'] = time.time() - start

    # - SCVI is trained unsupervised on all data to learn a latent space
    # - SCANVI is then trained in a semi-supervised fashion using the reference cell type labels
    # - Returns the combined AnnData object (`data`) and trained models (`vae`, `lvae`)

    start = time.time()
    data, vae, lvae = run_multi_model(gex_shared, query, ref)
    timing_summary['Model Training'] = time.time() - start

    # ----------------------------------------------
    # VAE – SCVI model (Single-cell Variational Inference)
    # ----------------------------------------------
    # SCVI is an unsupervised model that learns a low-dimensional latent representation
    # of the data by modeling the gene expression counts using a negative binomial
    # distribution with batch correction.
    # This model does not use cell type annotations (unsupervised).

    print(vae)

    # ----------------------------------------------
    # lvae – SCANVI model (Semi-supervised SCVI)
    # ----------------------------------------------
    # SCANVI is an extension of SCVI that incorporates label information (semi-supervised).
    # It uses the latent space learned by the VAE and propagates the known cell type labels
    # from the reference dataset to the unlabeled query data (GEX and ATAC).
    # This allows the model to predict the cell type labels for unlabeled cells.

    print(lvae)

    # AnnData object
    print(data)

    # Once the SCANVI model (lvae) is trained, we can use it to predict the cell types for
    # the unlabeled cells in the query dataset (ATAC and GEX). This function uses the
    # learned latent representation (from the VAE and SCANVI) to assign cell type labels
    # to these cells.
    #
    # Specifically, the model predicts the 'celltype_scanvi' labels for the query data
    # based on the reference annotations, and these predictions are stored in the `obs`
    # attribute of the AnnData object. The function also performs dimensionality reduction
    # (UMAP) based on the learned latent representation (`X_scANVI`) for visualization.
    #
    # The resulting `adata` object now includes predicted cell type labels and UMAP embeddings
    # that can be used for downstream analysis and visualization.
    #
    # The function outputs a UMAP plot and the modified AnnData object with the new labels.

    # Call the transfer_labels function to predict and visualize cell types in the data

    start = time.time()
    plot, data = transfer_labels(data, lvae)
    timing_summary['Label Transfer'] = time.time() - start

    # Write output
    data.write("atac_gex_scanvi")

    # ----------------------------------------------
    # Subset Multiome Data for RNA and ATAC
    # ----------------------------------------------
    # After the cell type annotations have been transferred to the query data, we now
    # need to separate the Multiome dataset (which contains both RNA and ATAC data) into
    # individual subsets for further analysis.
    #
    # This function takes the combined AnnData object (which contains both RNA and ATAC
    # modalities) and subsets it into two distinct datasets: one for the RNA data (`rna_unannotated`)
    # and another for the ATAC data (`atac_unannotated`). The subsetting is based on the "modality"
    # field, which differentiates between the two types of data.
    #
    # The following steps are performed:
    # - The `final_annotation` (predicted by SCANVI) is copied to each modality (RNA and ATAC).
    # - Common barcodes (cells) between RNA and ATAC datasets are ensured.
    # - Metadata related to biosample, donor, disease, organ, etc., is added to the `obs` attribute
    #   of both datasets.
    # - Raw counts are saved in the `.raw` attribute for future reference.
    #
    # The resulting two subsets are stored in separate AnnData objects, `rna_unannotated`
    # and `atac_unannotated`, and saved to disk for further use in downstream analyses.

    # Read h5ad output
    df = ad.read_h5ad("atac_gex_scanvi")
    # Print dataframe
    print(df)
    rna_unannotated, atac_unannotated = subset_multiome_data(df)
    # Check that we have raw counts
    is_integer = np.all(atac_unannotated.X == atac_unannotated.X.astype(int))
    # Check for barcodes between atac and gex datasets that do not share the same cell annotations, even though they come from the same cells
    mismatched = check_cell_matches(rna_unannotated, atac_unannotated)
    len(mismatched)
    # Cell types of unannotated RNA anndata
    rna_unannotated.obs.celltype
    # Cell type of unannotated ATAC anndata
    atac_unannotated.obs.celltype
    # Cell types of unannotated RNA anndata
    rna_unannotated.obs["celltype"]

    # ----------------------------------------------
    # Perform Label Transfer and Save Predictions
    # ----------------------------------------------
    # In this step, we transfer the predicted cell type labels from the SCANVI model
    # (`C_scANVI`) to the original `atac_shared` and `gex_shared` datasets.
    # The predictions are retrieved from the `data.obs` object, where cell type
    # labels are stored after training the SCANVI model.
    #
    # The labels are then saved back to the `final_annotation` field in both the
    # ATAC and GEX datasets.
    #
    # This step ensures that the cell type predictions from the reference
    # (annotated) dataset are applied to the query (unannotated) dataset.
    # The results will provide labeled annotations for both the ATAC and GEX
    # data matrices, allowing for further downstream analyses.
    # Save the annotations to gex and atac matrices
    atac_shared.obs['final_annotation'] = data.obs.loc[atac_shared.obs_names + '_atac_unannotated'][
        'C_scANVI'].to_numpy()
    gex_shared.obs['final_annotation'] = data.obs.loc[gex_shared.obs_names + '_rna_unannotated']['C_scANVI'].to_numpy()
    # Write the annotated matrices
    atac_shared.write("atac_annotated_matrix.h5ad")
    gex_shared.write("gex_annotated_matrix.h5ad")
    # Compare predicted cell type labels with Leiden clusters
    # The predicted annotations align well with Leiden cluster structure.
    # However, due to fewer cells in the ATAC-seq dataset, resolution is limited—
    # making it difficult to distinguish similar subtypes like CD4+ and CD8+ T cells,
    # as well as other closely related cell types.
    # Save SCANVI predictions for Multiome data using shared gene features
    # (Note: previous predictions may not have used a shared gene space across modalities)
    data.write("SCANVI_predictions.h5ad")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('-g', '--gex-file', required=True, help="Input gene expression AnnData file")
    parser.add_argument('-a', '--atac-file', required=True, help="Input ATAC AnnData file")
    parser.add_argument('-r', '--ref-file', required=True, help="Input label reference AnnData file")
    parser.add_argument(
        '-l',
        '--localize',
        required=True,
        default=False,
        action="store_true",
        help="Localize input files and push outputs back to bucket"
    )
    parsed_args = parser.parse_args()
    if parsed_args.localize:
        gex, atac, ref = pull_all_files([parsed_args.gex_file, parsed_args.atac_file, parsed_args.ref_file])
    else:
        gex, atac, ref = parsed_args.gex_file, parsed_args.atac_file, parsed_args.ref_file
    main(gex, atac, ref)
    if parsed_args.localize:
        bucket_name = get_bucket_and_path(parsed_args.ref_file)[0]
        delocalize_file(bucket_name, "SCANVI_predictions.h5ad", "SCANVI_predictions.h5ad")
