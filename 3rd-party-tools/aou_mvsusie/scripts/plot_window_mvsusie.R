#!/usr/bin/env Rscript

script_file <- sub(
  "^--file=",
  "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[[1L]]
)
script_file <- gsub("~\\+~", " ", script_file)
script_dir <- dirname(normalizePath(script_file))
source(file.path(script_dir, "trans_window_cli.R"))
source(file.path(script_dir, "trans_window_logging.R"))

args <- parse_cli_args(
  option_list = list(
    optparse::make_option("--prepared", type = "character"),
    optparse::make_option("--fit", type = "character"),
    optparse::make_option("--png", type = "character"),
    optparse::make_option("--pdf", type = "character"),
    optparse::make_option("--plot-rds", type = "character")
  ),
  description = "Plot joint mvSuSiE credible-set effects by feature."
)

prepared <- readRDS(require_cli_arg(args, "prepared"))
bundle <- readRDS(require_cli_arg(args, "fit"))
fit <- bundle$fit
if (!inherits(fit, "mvsusie")) {
  stop("The fit bundle does not contain an mvSuSiE fit.", call. = FALSE)
}

required_variant_columns <- c("CHROM", "POS", "variant_id")
if (!all(required_variant_columns %in% names(prepared$variant_metadata))) {
  stop("Prepared variant metadata is incomplete for plotting.", call. = FALSE)
}
chromosomes <- unique(sub("^chr", "", prepared$variant_metadata$CHROM))
chromosome <- suppressWarnings(as.integer(chromosomes))
if (length(chromosome) != 1L || is.na(chromosome)) {
  stop("The plot requires one numeric chromosome.", call. = FALSE)
}
positions_mb <- as.numeric(prepared$variant_metadata$POS) / 1e6
markers <- as.character(prepared$variant_metadata$variant_id)

required_outcome_columns <- c("outcome_key", "phenotype_id", "modality")
if (!all(required_outcome_columns %in% names(prepared$phenotype_metadata))) {
  stop("Prepared phenotype metadata is incomplete for plotting.", call. = FALSE)
}
display_labels <- paste0(
  prepared$phenotype_metadata$modality,
  ": ",
  prepared$phenotype_metadata$phenotype_id
)
if (anyDuplicated(display_labels)) {
  display_labels <- paste0(display_labels, " [", prepared$phenotype_metadata$outcome_key, "]")
}

pipeline_log("Calling mvsusieR::mvsusie_plot with conditional effects.")
plot_result <- mvsusieR::mvsusie_plot(
  fit = fit,
  chr = chromosome,
  pos = positions_mb,
  markers = markers,
  outcomes = display_labels,
  lfsr_cutoff = 0.05,
  sentinel_only = FALSE,
  add_cs = TRUE,
  conditional_effect = TRUE,
  sort_by_cs = TRUE
)
if (is.null(plot_result$effect_plot)) {
  stop("mvsusie_plot did not return an effect plot.", call. = FALSE)
}

png_path <- require_cli_arg(args, "png")
pdf_path <- require_cli_arg(args, "pdf")
plot_rds_path <- require_cli_arg(args, "plot_rds")
for (path in c(png_path, pdf_path, plot_rds_path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
}
ggplot2::ggsave(
  png_path,
  plot_result$effect_plot,
  width = 15,
  height = 10,
  dpi = 300
)
ggplot2::ggsave(
  pdf_path,
  plot_result$effect_plot,
  width = 15,
  height = 10
)
saveRDS(plot_result, plot_rds_path)
required_outputs <- c(png_path, pdf_path, plot_rds_path)
if (any(!file.exists(required_outputs)) || any(file.info(required_outputs)$size == 0)) {
  stop("The mvSuSiE plot outputs are incomplete.", call. = FALSE)
}
pipeline_log("The mvSuSiE credible-set-by-feature plots were saved.")
