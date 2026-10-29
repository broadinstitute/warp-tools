# 1.5.0
2026-10-29 (Date of Last Commit)

* Added scanvi_sample.py, the WARP entry point that labels one sample end to end, and validate_inputs.py, which checks every sample's files and ids before any sample runs
* Moved the WARP pipeline's per-sample preprocessing into preprocess_query.py (formerly inline in the WDL's PreprocessFilter task), unchanged
* Added scp_output.py: annotated matrices hold counts in .raw, log2(CPM + 1) in .X (units in bits, the Allen convention), a 2-column X_umap fitted on the sample's own latent and the 30-dimensional X_scANVI latent, written with gzip
* Added query mapping: a supplied model that does not know a sample's batch maps the sample onto it with scArches (load_query_data plus a short training on the sample); a supplied model is used without the reference, and the predictions file then holds the sample only
* Added run provenance to every output's uns: scanvi_input_id, scanvi_input_file, scanvi_model (the reference the model was trained on, from the model's reference.txt; reference_id.py) and scanvi_versions (key library versions)
* Removed main() and its helpers from multiome_label_transfer.py, gcs_utils.py and google-cloud-storage, the unused GENCODE v41 GFF3, and the scratch model saves, plots and data copies the training functions wrote
* Removed the scipy<1.14 pin (added 2026-09-04 for an anndata backed-CSR bug that no longer reproduces). The image now resolves anndata 0.13.4 and scipy 1.18.1, with scvi-tools 1.5.1, snapatac2 2.10.0, torch 2.15.0, lightning 2.6.6, scanpy 1.12.4, numpy 2.5.4 and pandas 2.3.3. Outputs are identical to the pinned build's on the WARP Plumbing data, and anndata 0.10.8 reads them
* Added /usr/local/versions.txt (the build's pip freeze), printed by the CI build so every build's versions are on record
* Added test_scp_and_query_mapping.py to CI, beside the seed test
* Published as us.gcr.io/broad-gotc-prod/scvi-scanvi@sha256:2f32323adcb5dafcca4ce4326316eb0a23dc9a1701ce96fd13b175e10907ffba, pinned by WARP scANVI 3.0.0

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
