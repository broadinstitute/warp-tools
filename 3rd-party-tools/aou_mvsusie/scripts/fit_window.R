#!/usr/bin/env Rscript

script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
script_file <- gsub("~\\+~", " ", script_file)
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "trans_window_io.R"))
source(file.path(script_dir, "trans_window_model.R"))
source(file.path(script_dir, "trans_window_prior.R"))
source(file.path(script_dir, "trans_window_cli.R"))

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--prepared", type = "character"),
    optparse::make_option("--window-id", type = "character"),
    optparse::make_option("--start-L", type = "integer", default = 10L),
    optparse::make_option("--step-L", type = "integer", default = 5L),
    optparse::make_option("--max-L", type = "integer", default = 40L),
    optparse::make_option(
      "--greedy-lbf-cutoff",
      type = "double",
      default = 1
    ),
    optparse::make_option("--max-iter", type = "integer", default = 100L),
    optparse::make_option("--tol", type = "double", default = 1e-4),
    optparse::make_option("--coverage", type = "double", default = 0.95),
    optparse::make_option("--min-abs-corr", type = "double", default = 0.5),
    optparse::make_option("--n-thread", type = "integer", default = 1L),
    optparse::make_option("--mashr-n-pca", type = "integer", default = 5L),
    optparse::make_option("--mashr-seed", type = "integer", default = NULL),
    optparse::make_option("--mashr-strong-lfsr", type = "double", default = 0.05),
    optparse::make_option("--mashr-output", type = "character"),
    optparse::make_option("--greedy-history-output", type = "character"),
    optparse::make_option("--covariate-provenance-output", type = "character"),
    optparse::make_option("--output", type = "character")
  ),
  description = "Fit mvSusieR for one prepared trans window."
)
pipeline_log("Reading prepared window data.")
prepared <- readRDS(require_cli_arg(args, "prepared"))
prepared <- validate_prepared_window(
  prepared,
  require_cli_arg(args, "window_id")
)
pipeline_log(sprintf(
  "Prepared data loaded: %d samples, %d variants, and %d outcomes.",
  nrow(prepared$X), ncol(prepared$X), ncol(prepared$Y)
))
mashr_seed <- optional_cli_arg(args, "mashr_seed")
if (!is.null(mashr_seed)) mashr_seed <- as_cli_integer(args, "mashr_seed", 0L)
config <- make_model_config(
  start_L = as_cli_integer(args, "start_L", 10L),
  step_L = as_cli_integer(args, "step_L", 5L),
  max_L = as_cli_integer(args, "max_L", 40L),
  greedy_lbf_cutoff = as_cli_numeric(args, "greedy_lbf_cutoff", 1),
  max_iter = as_cli_integer(args, "max_iter", 100L),
  tol = as_cli_numeric(args, "tol", 1e-4),
  coverage = as_cli_numeric(args, "coverage", 0.95),
  min_abs_corr = as_cli_numeric(args, "min_abs_corr", 0.5),
  n_thread = as_cli_integer(args, "n_thread", 1L),
  mashr_n_pca = as_cli_integer(args, "mashr_n_pca", 5L),
  mashr_seed = mashr_seed,
  mashr_strong_lfsr = as_cli_numeric(args, "mashr_strong_lfsr", 0.05)
)
result <- fit_window_mvsusie(prepared, config)
bundle <- list(
  fit = result$fit,
  metadata = result$metadata,
  window = prepared$window,
  variant_metadata = prepared$variant_metadata,
  phenotype_metadata = prepared$phenotype_metadata,
  samples = prepared$samples,
  covariate_rank = prepared$covariate_rank,
  phenotype_covariate_rank = prepared$phenotype_covariate_rank,
  qc = prepared$qc
)
pipeline_log("Saving the mvSuSiE fit bundle.")
save_rds_checked(bundle, require_cli_arg(args, "output"))
pipeline_log("Saving the mashr training bundle.")
save_rds_checked(
  result$mashr_training,
  require_cli_arg(args, "mashr_output")
)
pipeline_log("Saving the greedy L history.")
data.table::fwrite(
  result$greedy_history,
  require_cli_arg(args, "greedy_history_output"),
  sep = "\t"
)
provenance_output <- require_cli_arg(args, "covariate_provenance_output")
pipeline_log("Writing covariate provenance from the prepared window.")
data.table::fwrite(
  prepared$covariate_provenance,
  provenance_output,
  sep = "\t",
  quote = FALSE
)
if (!file.exists(provenance_output) || file.info(provenance_output)$size == 0) {
  stop("Failed to write covariate provenance: ", provenance_output, call. = FALSE)
}
pipeline_log("The mvSuSiE fit bundle was saved.")
