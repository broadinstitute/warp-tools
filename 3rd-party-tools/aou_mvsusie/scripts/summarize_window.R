#!/usr/bin/env Rscript

script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
script_file <- gsub("~\\+~", " ", script_file)
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "trans_window_model.R"))
source(file.path(script_dir, "trans_window_cli.R"))

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--prepared", type = "character"),
    optparse::make_option("--fit", type = "character"),
    optparse::make_option("--output-dir", type = "character")
  ),
  description = "Write per-window mvSusie posterior summary tables."
)
prepared <- readRDS(require_cli_arg(args, "prepared"))
bundle <- readRDS(require_cli_arg(args, "fit"))
output_dir <- require_cli_arg(args, "output_dir")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

window_id <- as.character(prepared$qc$window_id)
config <- bundle$metadata$config

pip <- extract_variant_pips(bundle$fit, prepared)
pip[, window_id := window_id]
pip <- merge(
  pip,
  data.table::copy(prepared$variant_metadata),
  by = "variant_id",
  all.x = TRUE,
  sort = FALSE
)
data.table::setcolorder(
  pip,
  c("window_id", "variant_id", "CHROM", "POS", "REF", "ALT", "pip")
)

credible_set_tables <- extract_credible_set_tables(bundle$fit, prepared, config)
credible_sets <- credible_set_tables$summary
credible_sets[, window_id := window_id]
data.table::setcolorder(
  credible_sets,
  c(
    "window_id", "component", "credible_set_size", "sentinel_variant_id",
    "sentinel_alpha", "coverage", "purity_min", "purity_mean"
  )
)

credible_set_members <- credible_set_tables$members
credible_set_members[, window_id := window_id]
data.table::setcolorder(
  credible_set_members,
  c("window_id", "component", "variant_id", "alpha", "pip", "is_sentinel")
)

component_feature_support <- extract_component_feature_support(bundle$fit, prepared)
component_feature_support[, window_id := window_id]
data.table::setcolorder(
  component_feature_support,
  c(
    "window_id", "component", "outcome_key", "modality", "phenotype_id",
    "single_effect_lfsr", "outcome_lbf"
  )
)

qc <- data.table::as.data.table(prepared$qc)
qc[, `:=`(
  converged = isTRUE(bundle$metadata$converged),
  niter = bundle$metadata$niter,
  n_credible_sets = length(bundle$fit$sets$cs),
  mvsusieR_version = bundle$metadata$mvsusieR_version,
  prior = bundle$metadata$prior,
  residual_variance_mode = bundle$metadata$residual_variance_mode,
  start_L = config$start_L,
  step_L = config$step_L,
  max_L = config$max_L,
  greedy_lbf_cutoff = config$greedy_lbf_cutoff,
  L_final = bundle$metadata$L_final
)]

data.table::fwrite(
  pip,
  file.path(output_dir, "variant_pip.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(
  credible_sets,
  file.path(output_dir, "credible_sets.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(
  credible_set_members,
  file.path(output_dir, "credible_set_members.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(
  component_feature_support,
  file.path(output_dir, "component_feature_support.tsv.gz"),
  sep = "\t",
  compress = "gzip"
)
data.table::fwrite(qc, file.path(output_dir, "window_qc.tsv"), sep = "\t")
