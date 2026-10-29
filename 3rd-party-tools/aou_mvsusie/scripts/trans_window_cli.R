parse_cli_args <- function(
    option_list = list(),
    args = commandArgs(trailingOnly = TRUE),
    description = "") {
  if (!requireNamespace("optparse", quietly = TRUE)) {
    stop("The optparse package is required for the command-line interface.", call. = FALSE)
  }
  parser <- optparse::OptionParser(
    option_list = option_list,
    description = description
  )
  parsed <- optparse::parse_args(
    parser,
    args = args,
    positional_arguments = FALSE,
    convert_hyphens_to_underscores = TRUE
  )
  parsed
}

require_cli_arg <- function(args, name) {
  value <- args[[name]]
  missing <- is.null(value) || !length(value) || identical(value, TRUE)
  if (is.character(value)) missing <- missing || !nzchar(value[[1L]])
  if (missing) {
    stop("Missing required command-line option: --", gsub("_", "-", name), call. = FALSE)
  }
  value
}

optional_cli_arg <- function(args, name, default = NULL) {
  value <- args[[name]]
  if (is.null(value) || identical(value, TRUE)) default else value
}

split_cli_paths <- function(value) {
  if (is.null(value) || !nzchar(value)) return(character())
  trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
}

as_cli_integer <- function(args, name, default) {
  value <- optional_cli_arg(args, name, as.character(default))
  out <- suppressWarnings(as.integer(value))
  if (is.na(out)) stop("--", name, " must be an integer.", call. = FALSE)
  out
}

as_cli_numeric <- function(args, name, default) {
  value <- optional_cli_arg(args, name, as.character(default))
  out <- suppressWarnings(as.numeric(value))
  if (is.na(out)) stop("--", name, " must be numeric.", call. = FALSE)
  out
}

save_rds_checked <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, path)
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write RDS output: ", path, call. = FALSE)
  }
}
