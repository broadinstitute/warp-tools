suppressPackageStartupMessages(library(data.table))

normalize_sample_ids <- function(ids) {
  ids <- trimws(as.character(ids))
  sub("^X(?=[0-9])", "", ids, perl = TRUE)
}

required_joint_modalities <- function() {
  c("expression", "splicing", "protein")
}

require_columns <- function(dt, required, label) {
  missing <- setdiff(required, names(dt))
  if (length(missing)) {
    stop(
      label, " is missing required columns: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

read_windows_manifest <- function(path) {
  dt <- fread(path, check.names = FALSE)
  require_columns(
    dt,
    c("window_id", "chrom", "start", "end", "dosage_file"),
    "Window manifest"
  )
  if (anyDuplicated(dt$window_id)) {
    stop("Window manifest contains duplicate window_id values.", call. = FALSE)
  }
  dt[, start := as.integer(start)]
  dt[, end := as.integer(end)]
  if (anyNA(dt$start) || anyNA(dt$end) || any(dt$start < 0L) || any(dt$end <= dt$start)) {
    stop("Window coordinates must be valid 0-based half-open intervals.", call. = FALSE)
  }
  dt[]
}

read_window_phenotypes_manifest <- function(path) {
  dt <- fread(path, check.names = FALSE)
  require_columns(
    dt,
    c("window_id", "outcome_key", "phenotype_id", "modality", "phenotype_file"),
    "Window phenotype manifest"
  )
  if (!nrow(dt)) {
    stop("Window phenotype manifest must contain at least one outcome.", call. = FALSE)
  }
  if (anyDuplicated(dt[, paste(window_id, outcome_key, sep = "\r")])) {
    stop(
      "Window phenotype manifest contains duplicate window_id/outcome_key pairs.",
      call. = FALSE
    )
  }
  expected_outcome_keys <- paste(dt$modality, dt$phenotype_id, sep = "::")
  if (any(dt$outcome_key != expected_outcome_keys)) {
    stop(
      "Each outcome_key must equal modality::phenotype_id.",
      call. = FALSE
    )
  }
  unsupported <- setdiff(
    unique(as.character(dt$modality)),
    required_joint_modalities()
  )
  if (length(unsupported)) {
    stop(
      "Window phenotype manifest contains an unsupported modality: ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }
  dt[]
}

numeric_matrix <- function(values, label) {
  out <- suppressWarnings(matrix(as.numeric(as.matrix(values)), nrow = nrow(values)))
  if (any(!is.finite(out) & !is.na(out))) {
    stop(label, " contains non-numeric values.", call. = FALSE)
  }
  out
}

validate_prepared_window <- function(prepared, expected_window_id) {
  required <- c(
    "window", "X", "Y", "variant_metadata", "phenotype_metadata",
    "samples", "covariate_provenance", "covariate_rank",
    "phenotype_covariate_rank", "qc"
  )
  missing <- setdiff(required, names(prepared))
  if (length(missing)) {
    stop(
      "Prepared window is missing required fields: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  actual_window_id <- as.character(prepared$window$window_id)
  if (
    length(actual_window_id) != 1L || is.na(actual_window_id) ||
    actual_window_id != expected_window_id
  ) {
    stop(
      "Prepared window ID does not match requested window: ",
      expected_window_id,
      call. = FALSE
    )
  }
  if (
    !is.matrix(prepared$X) || !is.matrix(prepared$Y) ||
    nrow(prepared$X) != nrow(prepared$Y) ||
    nrow(prepared$X) != length(prepared$samples)
  ) {
    stop("Prepared X, Y, and samples must be aligned.", call. = FALSE)
  }
  if (
    ncol(prepared$X) != nrow(prepared$variant_metadata) ||
    ncol(prepared$Y) != nrow(prepared$phenotype_metadata)
  ) {
    stop("Prepared matrices and metadata must have matching sizes.", call. = FALSE)
  }
  if (!is.data.frame(prepared$covariate_provenance)) {
    stop("Prepared covariate provenance must be a data frame.", call. = FALSE)
  }
  prepared
}

read_wide_dosage <- function(path) {
  dt <- fread(path, check.names = FALSE)
  required <- c("CHROM", "POS", "REF", "ALT")
  require_columns(dt, required, "Dosage file")
  sample_columns <- setdiff(names(dt), required)
  if (!length(sample_columns)) {
    stop("Dosage file has no sample columns.", call. = FALSE)
  }
  sample_ids <- normalize_sample_ids(sample_columns)
  if (anyDuplicated(sample_ids)) {
    stop("Dosage file contains duplicate sample IDs.", call. = FALSE)
  }
  X <- numeric_matrix(dt[, ..sample_columns], "Dosage values")
  X <- t(X)
  rownames(X) <- sample_ids
  variant_ids <- paste0(dt$CHROM, ":", dt$POS, "_", dt$REF, "_", dt$ALT)
  if (anyDuplicated(variant_ids)) {
    stop("Dosage file contains duplicate variant IDs.", call. = FALSE)
  }
  colnames(X) <- variant_ids
  list(
    metadata = dt[, ..required],
    X = X,
    sample_ids = sample_ids,
    variant_ids = variant_ids
  )
}

phenotype_layout <- function(modality) {
  if (modality %in% required_joint_modalities()) {
    return(list(id_column = 4L, metadata_columns = 1:4))
  }
  stop("Unsupported phenotype modality: ", modality, call. = FALSE)
}

read_phenotype_rows <- function(path, modality, phenotype_ids) {
  dt <- fread(path, check.names = FALSE)
  layout <- phenotype_layout(modality)
  if (ncol(dt) <= max(layout$metadata_columns)) {
    stop("Phenotype file has no sample columns: ", path, call. = FALSE)
  }
  id_values <- as.character(dt[[layout$id_column]])
  if (anyDuplicated(id_values)) {
    stop("Phenotype file contains duplicate feature IDs: ", path, call. = FALSE)
  }
  sample_columns <- setdiff(names(dt), names(dt)[layout$metadata_columns])
  sample_ids <- normalize_sample_ids(sample_columns)
  if (anyDuplicated(sample_ids)) {
    stop("Phenotype file contains duplicate sample IDs: ", path, call. = FALSE)
  }
  phenotype_ids <- as.character(phenotype_ids)
  missing_ids <- setdiff(phenotype_ids, id_values)
  if (length(missing_ids)) {
    stop(
      "Missing phenotype IDs in ", path, ": ",
      paste(missing_ids, collapse = ", "),
      call. = FALSE
    )
  }
  row_index <- match(phenotype_ids, id_values)
  values <- numeric_matrix(dt[row_index, ..sample_columns], "Phenotype values")
  Y <- t(values)
  rownames(Y) <- sample_ids
  colnames(Y) <- phenotype_ids
  metadata_columns <- names(dt)[layout$metadata_columns]
  metadata <- dt[row_index, ..metadata_columns]
  metadata[, modality := modality]
  list(
    Y = Y,
    metadata = metadata,
    sample_ids = sample_ids,
    phenotype_ids = phenotype_ids,
    modality = modality,
    source_file = path
  )
}

resolve_file_reference <- function(reference, files) {
  files <- as.character(files)
  direct <- files[files == reference]
  if (length(direct) == 1L) return(direct)
  matches <- files[basename(files) == basename(reference)]
  if (length(matches) != 1L) {
    stop(
      "Could not resolve file reference ", reference,
      " among supplied files.",
      call. = FALSE
    )
  }
  matches
}

read_window_phenotypes <- function(window_id, phenotype_manifest, phenotype_files) {
  rows <- phenotype_manifest[phenotype_manifest[["window_id"]] == window_id]
  if (!nrow(rows)) {
    stop("No phenotype rows found for window: ", window_id, call. = FALSE)
  }
  grouped_rows <- split(seq_len(nrow(rows)), rows$phenotype_file)
  parts <- lapply(grouped_rows, function(indices) {
    source_file <- resolve_file_reference(rows$phenotype_file[[indices[[1L]]]], phenotype_files)
    modalities <- as.character(rows$modality[indices])
    data <- read_phenotype_rows(
      source_file,
      modalities[[1L]],
      rows$outcome_key[indices]
    )
    data$metadata[, modality := modalities]
    list(indices = indices, data = data)
  })
  common_samples <- Reduce(intersect, lapply(parts, function(part) part$data$sample_ids))
  if (!length(common_samples)) {
    stop("No shared phenotype samples for window: ", window_id, call. = FALSE)
  }
  Y <- do.call(cbind, lapply(parts, function(part) {
    part$data$Y[common_samples, , drop = FALSE]
  }))
  column_order <- order(unlist(lapply(parts, `[[`, "indices")))
  Y <- Y[, column_order, drop = FALSE]
  metadata <- copy(rows)
  list(
    Y = Y,
    metadata = metadata,
    sample_ids = common_samples,
    phenotype_ids = as.character(rows$outcome_key),
    modalities = as.character(rows$modality),
    source_files = rows$phenotype_file
  )
}

read_covariate_file <- function(path) {
  dt <- fread(path, check.names = FALSE)
  if (
    identical(names(dt)[[1L]], "V1") && nrow(dt) > 0L &&
    identical(toupper(as.character(dt[[1L]][[1L]])), "ID")
  ) {
    header <- as.character(dt[1L])
    dt <- dt[-1L]
    setnames(dt, header)
  }
  if (ncol(dt) < 2L) stop("Covariate file has no sample columns: ", path, call. = FALSE)
  covariate_ids <- as.character(dt[[1L]])
  if (anyDuplicated(covariate_ids)) {
    stop("Covariate file contains duplicate covariate IDs: ", path, call. = FALSE)
  }
  sample_columns <- names(dt)[-1L]
  sample_ids <- normalize_sample_ids(sample_columns)
  if (anyDuplicated(sample_ids)) {
    stop("Covariate file contains duplicate sample IDs: ", path, call. = FALSE)
  }
  values <- numeric_matrix(dt[, -1L, with = FALSE], "Covariate values")
  values <- t(values)
  rownames(values) <- sample_ids
  colnames(values) <- covariate_ids
  values
}

read_joint_covariates <- function(
    expression_path,
    splicing_path,
    protein_path
) {
  paths <- c(
    expression = expression_path,
    splicing = splicing_path,
    protein = protein_path
  )
  if (any(!file.exists(paths))) {
    stop("Every joint covariate file must exist.", call. = FALSE)
  }
  result <- lapply(paths, read_covariate_file)
  shared_samples <- Reduce(intersect, lapply(result, rownames))
  if (!length(shared_samples)) {
    stop("Joint covariate files have no shared sample IDs.", call. = FALSE)
  }
  result
}

unique_covariate_columns <- function(matrices) {
  if (!length(matrices)) stop("At least one covariate matrix is required.", call. = FALSE)
  sample_ids <- rownames(matrices[[1L]])
  if (is.null(sample_ids)) stop("Covariate matrices must have sample IDs.", call. = FALSE)
  output <- matrix(numeric(), nrow = length(sample_ids), ncol = 0L)
  rownames(output) <- sample_ids
  seen <- character()
  for (matrix in matrices) {
    if (!identical(rownames(matrix), sample_ids)) {
      stop("Covariate matrices do not contain the same sample IDs.", call. = FALSE)
    }
    duplicate_ids <- intersect(seen, colnames(matrix))
    for (covariate_id in duplicate_ids) {
      if (!isTRUE(all.equal(output[, covariate_id], matrix[, covariate_id]))) {
        stop(
          "Covariate ID has conflicting values across files: ",
          covariate_id,
          call. = FALSE
        )
      }
    }
    new_ids <- setdiff(colnames(matrix), seen)
    if (length(new_ids)) {
      output <- cbind(output, matrix[, new_ids, drop = FALSE])
      seen <- c(seen, new_ids)
    }
  }
  output
}
