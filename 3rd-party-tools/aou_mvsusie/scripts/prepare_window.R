#!/usr/bin/env Rscript

script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
script_file <- gsub("~\\+~", " ", script_file)
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "trans_window_io.R"))
source(file.path(script_dir, "trans_window_logging.R"))
source(file.path(script_dir, "trans_window_preprocess.R"))
source(file.path(script_dir, "trans_window_cli.R"))

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--windows", type = "character"),
    optparse::make_option("--window-phenotypes", type = "character"),
    optparse::make_option("--window-id", type = "character"),
    optparse::make_option("--dosage", type = "character"),
    optparse::make_option("--phenotype-files", type = "character"),
    optparse::make_option("--expression-covariates", type = "character"),
    optparse::make_option("--splicing-covariates", type = "character"),
    optparse::make_option("--protein-covariates", type = "character"),
    optparse::make_option("--keep-samples", type = "character", default = NULL),
    optparse::make_option("--min-genotype-variance", type = "double", default = 1e-8),
    optparse::make_option("--min-phenotype-variance", type = "double", default = 1e-8),
    optparse::make_option("--covariate-provenance-output", type = "character"),
    optparse::make_option("--output", type = "character")
  ),
  description = "Prepare one trans-window mvSusie input bundle."
)
pipeline_log("Reading window and phenotype manifests.")
windows <- read_windows_manifest(require_cli_arg(args, "windows"))
phenotype_manifest <- read_window_phenotypes_manifest(
  require_cli_arg(args, "window_phenotypes")
)
window_id <- require_cli_arg(args, "window_id")
window <- windows[windows[["window_id"]] == window_id]
if (nrow(window) != 1L) {
  stop("Expected exactly one window row for window_id: ", window_id, call. = FALSE)
}

phenotype_files <- split_cli_paths(require_cli_arg(args, "phenotype_files"))
pipeline_log("Reading genotype data.")
dosage <- read_wide_dosage(require_cli_arg(args, "dosage"))
pipeline_log(sprintf(
  "Genotype data loaded: %d samples and %d variants.",
  nrow(dosage$X), ncol(dosage$X)
))
pipeline_log("Reading selected phenotype data.")
phenotype_data <- read_window_phenotypes(window_id, phenotype_manifest, phenotype_files)
pipeline_log(sprintf(
  "Phenotype data loaded: %d samples and %d outcomes.",
  nrow(phenotype_data$Y), ncol(phenotype_data$Y)
))
pipeline_log("Reading modality-specific covariates.")
covariates_by_modality <- read_joint_covariates(
  expression_path = require_cli_arg(args, "expression_covariates"),
  splicing_path = require_cli_arg(args, "splicing_covariates"),
  protein_path = require_cli_arg(args, "protein_covariates")
)
for (modality in names(covariates_by_modality)) {
  pipeline_log(sprintf(
    "%s covariates loaded: %d samples and %d columns.",
    modality,
    nrow(covariates_by_modality[[modality]]),
    ncol(covariates_by_modality[[modality]])
  ))
}

pipeline_log("Residualizing genotype and phenotype matrices.")
prepared <- prepare_joint_window_data(
  window = window,
  phenotype_data = phenotype_data,
  dosage = dosage,
  covariates_by_modality = covariates_by_modality,
  keep_samples = optional_cli_arg(args, "keep_samples"),
  min_genotype_variance = as_cli_numeric(args, "min_genotype_variance", 1e-8),
  min_phenotype_variance = as_cli_numeric(args, "min_phenotype_variance", 1e-8)
)
prepared$input_checksums <- input_file_checksums(c(
  windows = require_cli_arg(args, "windows"),
  window_phenotypes = require_cli_arg(args, "window_phenotypes"),
  dosage = require_cli_arg(args, "dosage"),
  stats::setNames(
    phenotype_files,
    paste0("phenotype_file_", seq_along(phenotype_files))
  ),
  expression_covariates = require_cli_arg(args, "expression_covariates"),
  splicing_covariates = require_cli_arg(args, "splicing_covariates"),
  protein_covariates = require_cli_arg(args, "protein_covariates")
))
pipeline_log(sprintf(
  "Preprocessing complete: %d samples, %d variants, and %d outcomes retained.",
  nrow(prepared$X), ncol(prepared$X), ncol(prepared$Y)
))

pipeline_log("Saving prepared window data.")
save_rds_checked(prepared, require_cli_arg(args, "output"))
pipeline_log("Writing covariate provenance.")
write_covariate_provenance(
  prepared$covariate_provenance,
  require_cli_arg(args, "covariate_provenance_output")
)
pipeline_log("Prepared window data saved.")
