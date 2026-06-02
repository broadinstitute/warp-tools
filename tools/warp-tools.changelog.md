# 2.7.1
2026-06-01 (Date of Last Commit)

* Added `--input_id_name` argument to `create_h5ad_optimus.py` and `create_snrna_optimus_full_h5ad.py`; this optional argument (default `"input_id"`) controls the metadata key name under which the `input_id` value is stored in h5ad `obs` columns and `uns` global attributes

# 2.7.0
2026-01-26 (Date of Last Commit)

* Added script test_metrics_zeros_multiome.py which tests metrics calculation on multiome data for zero values
# 2.6.1
2025-02-12 (Date of Last Commit)

* Updated the sctools repository to use concat instead of append, since append is now deprecated