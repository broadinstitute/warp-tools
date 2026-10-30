# 1.7.0
2026-10-30 (Date of Last Commit)

* Added gene matching by Ensembl ID (a SCORCH sample has 99.5 % of AIT33's genes by Ensembl ID, 93.3 % by symbol): when a supplied model's genes are Ensembl IDs, the sample's genes are looked up in its Ensembl ID column (Optimus's ensembl_ids) before matching; the annotated outputs keep the sample's gene names. Every supplied-model run prints how many of the model's genes the sample has
* Published as us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:59206dd15b3fb8170b63c0da429b35ba597847ee30e4cf520b513ed9bebad22e, pinned by WARP scANVI 3.0.0

# 1.6.0
2026-10-30 (Date of Last Commit)

* Renamed obsm['X_scANVI_<level>'] to obsm['X_scVI_<level>']: the 30-dimensional latent is the scVI encoder's z (as trained with that level's SCANVI model), not a SCANVI-specific space. Values are unchanged
* Published as us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:938baa80261c9d2122838403d36e000cf914035705f4cbddc6eed34595289386

# 1.5.0
2026-10-29 (Date of Last Commit)

* Added scanvi_sample.py, the WARP entry point that labels one sample end to end, and validate_inputs.py, which checks every sample's files and ids, and that no two outputs share a name, before any sample runs
* Added label levels: a sample is labelled at every level of the model, one SCANVI model per level (train_scvi once, then train_scanvi per level when training). Each level gives obs['<level>'], obsm['X_scANVI_<level>'] and obsm['X_umap_<level>'] (and <level>_max_probability); X_umap is a copy of the first level's and uns['scanvi_levels'] lists them. run_multi_model, run_gex_only_model, _concat_and_train_models, transfer_labels and finalize_output are replaced, and the placeholder final_annotation column is gone
* Added model_archive.py (replacing reference_id.py): one model archive holds a reference's scVI model and its SCANVI model per level, with reference.txt and levels.txt; label_transfer() reads and writes it, and refuses a single-model archive
* Changed the outputs to the paths given by --gex-out, --atac-out and --model-out (the model only when training) instead of names built from --input-id, and removed the predictions file
* Moved the WARP pipeline's per-sample preprocessing into preprocess_query.py (formerly inline in the WDL's PreprocessFilter task), unchanged
* Added scp_output.py: annotated matrices hold counts in .raw, log2(CPM + 1) in .X (units in bits, the Allen convention), and per level a 2-column UMAP fitted on the sample's own latent and the 30-dimensional scANVI latent, written with gzip
* Added query mapping: a supplied model that does not know a sample's batch maps the sample onto it with scArches (load_query_data plus a short training on the sample); a supplied model is used without the reference
* Added run provenance to every output's uns: scanvi_input_id, scanvi_input_file, scanvi_model (the reference the model was trained on, from the model archive's reference.txt), scanvi_levels and scanvi_versions (key library versions)
* Removed main() and its helpers from multiome_label_transfer.py, gcs_utils.py and google-cloud-storage, the unused GENCODE v41 GFF3, and the scratch model saves, plots and data copies the training functions wrote
* Removed the scipy<1.14 pin (added 2026-09-04 for an anndata backed-CSR bug that no longer reproduces). The image now resolves anndata 0.13.4 and scipy 1.18.1, with scvi-tools 1.5.1, snapatac2 2.10.0, torch 2.15.0, lightning 2.6.6, scanpy 1.12.5, numpy 2.5.4 and pandas 2.3.3. Outputs are identical to the pinned build's on the WARP Plumbing data, and anndata 0.10.8 reads them
* Added /usr/local/versions.txt (the build's pip freeze), printed by the CI build so every build's versions are on record
* Added test_scp_and_query_mapping.py to CI, beside the seed test
* Updated make_reference_only_model.py to train the way warp-scanvi-training does since 2026-10-29: negative binomial, dispersion gene, scvi-tools' scArches settings (covariates encoded, layer norm, no batch norm, dropout 0.2) and scVI set up without labels, then a SCANVI per --levels entry, into a model archive; default 200 epochs, because those settings converge slowly on a small reference
* Published as us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:eaadce983e21a8cac8461d99bc20489b887451f0e24f364a7f9a1013b972d811

# 1.4.0
2026-10-22 (Date of Last Commit)

* Removed the twelve placeholder sample-metadata columns that finalize_output stamped onto the predictions file
* Built with scvi-tools 1.5.1, snapatac2 2.10.0, anndata 0.12.19 and scipy below 1.14

# 1.3.0
2026-10-09 (Date of Last Commit)

* Added --seed for reproducible SCVI/SCANVI training and the test_seed_reproducibility.py CI check

# 1.2.0
2026-09-19 (Date of Last Commit)

* Updated to scvi-tools 1.5.1 and snapatac2 2.10.0 (from 1.2.0 and 2.7.0)
* Required exactly one model.pt in a supplied model archive and extracted it with Python's safe data filter

# 1.1.0
2026-06-27 (Date of Last Commit)

* Added run_gex_only_model for GEX-only label transfer

# 1.0.0
2025-08-26 (Date of Last Commit)

* Initial image: scvi-tools 1.2 and snapatac2 2.7 for multiome SCVI/SCANVI label transfer
