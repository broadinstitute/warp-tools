import numpy as np
import anndata as ad
import pandas as pd
from sklearn.neighbors import KNeighborsTransformer
import argparse
import scanpy as sc


def compute_doublet_scores(gex_h5ad, proportion_artificial=0.2):
    adata = ad.read_h5ad(gex_h5ad)
    adata.var_names_make_unique()
    k = np.int64(np.round(np.min([100, adata.shape[0] * 0.01])))
    n_doublets = np.int64(np.round(adata.shape[0] / (1 - proportion_artificial) - adata.shape[0]))
    real_cells_1 = np.random.choice(adata.obs_names, size=n_doublets, replace=True)
    real_cells_2 = np.random.choice(adata.obs_names, size=n_doublets, replace=True)
    doublet_X = adata[real_cells_1, :].X + adata[real_cells_2, :].X
    doublet_obs_names = ["X" + str(i) for i in range(n_doublets)]
    doublet_adata = ad.AnnData(X=doublet_X, obs=pd.DataFrame(index=doublet_obs_names), var=pd.DataFrame(index=adata.var_names))    
    adata = adata.concatenate(doublet_adata, index_unique=None)

    adata.obs["doublet_cell"] = adata.obs_names.isin(doublet_obs_names)
    adata.obs["doublet_cell"] = adata.obs["doublet_cell"].astype("category")
    adata.layers["UMIs"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=5000, flavor="seurat_v3", layer="UMIs")
        adata.layers["UMIs"]

    except:
        sc.pp.highly_variable_genes(adata, min_mean=1, min_disp=0.5)
        del adata.layers["UMIs"]

    adata_sub = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata_sub)
    sc.pp.pca(adata_sub)
    v = adata_sub.uns['pca']['variance']
    n_pcs = np.max(np.where(((v - np.mean(v)) / np.std(v)) > 0)[0])   
    knn = KNeighborsTransformer(
        n_neighbors=k,
        algorithm="kd_tree",
        n_jobs=-1,
    )

    knn = knn.fit(adata_sub.obsm["X_pca"][:, :n_pcs])
    dist, idx = knn.kneighbors()  
    knn_mapper = KNeighborsTransformer(
        n_neighbors=10,
        algorithm="kd_tree",
        n_jobs=-1,
    )

    knn_mapper = knn_mapper.fit(adata_sub[adata_sub.obs["doublet_cell"] == False, :].obsm["X_pca"][:, :n_pcs])
    dist1, _ = knn_mapper.kneighbors(adata_sub[adata_sub.obs["doublet_cell"] == True, :].obsm["X_pca"][:, :n_pcs])    
    dist_th = np.mean(dist1) + (1.64 * np.std(dist1))
    freq = (dist < dist_th) & (idx > adata[adata.obs["doublet_cell"] == False, :].shape[0])
    score1 = freq.mean(axis=1)
    score2 = freq[:, :np.int64(np.ceil(k/2))].mean(axis=1)
    adata.obs["doublet_score"] = [np.max([score1[i], score2[i]]) for i in range(adata.shape[0])]   
    end_csv=adata.obs.loc[~adata.obs_names.isin(doublet_obs_names), ["doublet_score"]]
    
    # Calculate the percentage of doublets with a doublet_score > 0.3
    num_doublets = end_csv[end_csv["doublet_score"] > 0.3].shape[0]
    total_cells = end_csv.shape[0]
    percent_doublets = num_doublets / total_cells * 100

    return end_csv, percent_doublets 


def process_gex_data(gex_h5ad, gex_nhash_id, library_csv, input_id, doublets):
    print("Reading Optimus h5ad:")
    gex_data = ad.read_h5ad(gex_h5ad)
    gex_data.uns['NHashID'] = gex_nhash_id
    gex_data.write(f"{input_id}.h5ad")

    print("Reading library metrics")
    library = pd.read_csv(library_csv, header=None)

    print("Calculating TSO frac")
    tso_reads = gex_data.obs.tso_reads.sum() / gex_data.obs.n_reads.sum()
    print("TSO reads:")
    print(tso_reads)
    dictionary = library.set_index(0)[1].to_dict()
    dictionary['frac_tso'] = tso_reads
    dictionary['percent_doublets'] = doublets
    new_dictionary = {"NHashID": [gex_nhash_id]}  # This line is fine, it already has a list
    # Update other scalar values to lists
    dictionary = {key: [value] for key, value in dictionary.items()}
    new_dictionary.update(dictionary)
    new_dictionary = pd.DataFrame(new_dictionary)
    new_dictionary.transpose().to_csv(f"{input_id}_{gex_nhash_id}_library_metrics.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process single-cell RNA-seq data and compute doublet scores.")
    parser.add_argument("--proportion_artificial", type=float, default=0.2, help="Proportion of artificial doublets to be generated (default is 0.2).")
    parser.add_argument("--gex_h5ad", type=str, required=True, help="Path to the GEX h5ad file.")
    parser.add_argument("--gex_nhash_id", type=str, required=True, help="NHashID for the GEX data.")
    parser.add_argument("--library_csv", type=str, required=True, help="Path to the library metrics CSV file.")
    parser.add_argument("--input_id", type=str, required=True, help="Input ID for output files.")

    args = parser.parse_args()

    # Load the AnnData object

    # Compute doublet scores
    doublet_scores, percent_doublets = compute_doublet_scores(args.gex_h5ad, proportion_artificial=args.proportion_artificial)
    print("Printing adata: ", doublet_scores)

    # Output the results
    output_path = args.gex_h5ad.replace(".h5ad", "_doublet_scores.csv")
    print("Output path: ", output_path)
    doublet_scores.to_csv(output_path)
    print(f"Doublet scores saved to {output_path}")
    print("Percent_doublets: ", percent_doublets)

    # Process GEX data
    process_gex_data(args.gex_h5ad, args.gex_nhash_id, args.library_csv, args.input_id, percent_doublets)
