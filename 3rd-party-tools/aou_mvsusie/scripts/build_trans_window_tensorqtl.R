#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

tensorqtl_columns <- c(
  "variant_id", "phenotype_id", "pval", "b", "b_se", "af"
)

require_columns <- function(data, required, label) {
  available <- if (is.character(data)) data else names(data)
  missing <- setdiff(required, available)

  if (length(missing) > 0L) {
    stop(
      label,
      " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

format_coordinate <- function(value) {
  format(value, scientific = FALSE, trim = TRUE)
}

parse_variant_locations <- function(associations, window_size_bp) {
  associations %>%
    mutate(
      chrom = str_extract(variant_id, "^[^:]+"),
      pos = as.integer(str_extract(variant_id, "(?<=:)[0-9]+"))
    ) %>%
    filter(!is.na(chrom), !is.na(pos), pos > 0L) %>%
    mutate(
      window_index = floor((pos - 1) / window_size_bp),
      start = window_index * window_size_bp,
      end = (window_index + 1) * window_size_bp,
      window_id = paste(
        chrom,
        format_coordinate(start),
        format_coordinate(end),
        sep = "_"
      )
    )
}

build_trans_window_tensorqtl_outputs <- function(
    trans_associations,
    output_dir,
    window_size_bp = 2e6,
    trans_p_threshold = 1e-8
) {
  require_columns(
    trans_associations,
    c(tensorqtl_columns, "modality"),
    "Trans TensorQTL associations"
  )

  if (length(window_size_bp) != 1L ||
      !is.finite(window_size_bp) || window_size_bp <= 0) {
    stop("window_size_bp must be a positive finite number.", call. = FALSE)
  }
  if (length(trans_p_threshold) != 1L ||
      !is.finite(trans_p_threshold) || trans_p_threshold <= 0) {
    stop("trans_p_threshold must be a positive finite number.", call. = FALSE)
  }

  trans <- trans_associations %>%
    mutate(.p_value = suppressWarnings(as.numeric(.data$pval))) %>%
    filter(
      !is.na(.data$.p_value),
      is.finite(.data$.p_value),
      .data$.p_value < trans_p_threshold
    ) %>%
    parse_variant_locations(window_size_bp)

  if (nrow(trans) == 0L) {
    stop("No trans associations pass the p-value threshold.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  associations <- trans %>%
    transmute(
      window_id = .data$window_id,
      chrom = .data$chrom,
      start = .data$start,
      end = .data$end,
      modality = as.character(.data$modality),
      molecular_trait_id = as.character(.data$phenotype_id),
      p_value = .data$.p_value
    ) %>%
    group_by(
      .data$window_id,
      .data$chrom,
      .data$start,
      .data$end,
      .data$modality,
      .data$molecular_trait_id
    ) %>%
    summarise(p_value = min(.data$p_value), .groups = "drop") %>%
    arrange(.data$chrom, .data$start, .data$modality, .data$p_value, .data$molecular_trait_id)

  associations_path <- file.path(
    output_dir,
    "trans_window_associations.tsv.gz"
  )
  write_tsv(associations, associations_path)

  qc <- associations %>%
    group_by(
      .data$window_id,
      .data$chrom,
      .data$start,
      .data$end
    ) %>%
    summarise(
      n_trait_mappings = dplyr::n(),
      n_expression = sum(.data$modality == "expression"),
      n_splicing = sum(.data$modality == "splicing"),
      n_protein = sum(.data$modality == "protein"),
      min_p_value = min(.data$p_value),
      .groups = "drop"
    ) %>%
    arrange(.data$chrom, .data$start)
  qc_path <- file.path(output_dir, "trans_window_qc.tsv.gz")
  write_tsv(qc, qc_path)

  list(
    associations = associations,
    associations_path = normalizePath(associations_path, mustWork = TRUE),
    qc = qc,
    qc_path = normalizePath(qc_path, mustWork = TRUE)
  )
}

read_trans_input_files <- function(
    paths,
    labels,
    trans_p_threshold = Inf
) {
  paths <- as.character(paths)
  labels <- as.character(labels)

  if (length(paths) == 0L || length(paths) != length(labels)) {
    stop(
      "The number of trans input files must match the number of labels.",
      call. = FALSE
    )
  }
  if (anyNA(paths) || any(!nzchar(paths))) {
    stop("Trans input file paths cannot be empty.", call. = FALSE)
  }
  if (anyNA(labels) || any(!nzchar(labels))) {
    stop("Trans input labels cannot be empty.", call. = FALSE)
  }
  if (length(trans_p_threshold) != 1L || is.na(trans_p_threshold) ||
      trans_p_threshold <= 0) {
    stop("trans_p_threshold must be a positive number.", call. = FALSE)
  }

  map2_dfr(paths, labels, function(path, modality) {
    message("Reading ", modality, " trans associations from ", path, ".")
    if (str_detect(path, regex("[.]parquet$", ignore_case = TRUE))) {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop(
          "The arrow R package is required to read Parquet trans associations.",
          call. = FALSE
        )
      }
      dataset <- arrow::open_dataset(path, format = "parquet")
      require_columns(dataset$schema$names, tensorqtl_columns, "Parquet input")
      message(
        "Applying the Parquet p-value filter before collecting ",
        modality,
        " rows."
      )
      data <- dataset %>%
        select(all_of(tensorqtl_columns)) %>%
        filter(.data$pval < trans_p_threshold) %>%
        collect()
    } else {
      data <- read_tsv(path, show_col_types = FALSE)
      if (is.finite(trans_p_threshold)) {
        data <- data %>%
          mutate(pval = suppressWarnings(as.numeric(.data$pval))) %>%
          filter(
            !is.na(.data$pval),
            is.finite(.data$pval),
            .data$pval < trans_p_threshold
          )
      }
    }
    message("Retained ", nrow(data), " ", modality, " trans rows.")
    data %>% mutate(modality = modality)
  })
}

main <- function() {
  source("scripts/trans_window_cli.R")
  args <- parse_cli_args(
    option_list = list(
      optparse::make_option("--trans-files", type = "character"),
      optparse::make_option("--trans-labels", type = "character"),
      optparse::make_option("--output-dir", type = "character"),
      optparse::make_option("--window-size-bp", type = "double", default = 2e6),
      optparse::make_option("--trans-p-threshold", type = "double", default = 1e-8)
    ),
    description = "Build trans TensorQTL files by genomic window."
  )

  trans_files <- split_cli_paths(require_cli_arg(args, "trans_files"))
  trans_labels <- split_cli_paths(require_cli_arg(args, "trans_labels"))

  trans_associations <- read_trans_input_files(
    paths = trans_files,
    labels = trans_labels,
    trans_p_threshold = as.numeric(args$trans_p_threshold)
  )

  result <- build_trans_window_tensorqtl_outputs(
    trans_associations = trans_associations,
    output_dir = require_cli_arg(args, "output_dir"),
    window_size_bp = as.numeric(args$window_size_bp),
    trans_p_threshold = as.numeric(args$trans_p_threshold)
  )

  message("Wrote ", nrow(result$associations), " trans-window molecular-trait mappings.")
  message("Associations: ", result$associations_path)
  message("Window QC: ", result$qc_path)
}

if (sys.nframe() == 0L) {
  main()
}
