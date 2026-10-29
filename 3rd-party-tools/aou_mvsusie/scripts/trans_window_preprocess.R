rank_int <- function(x) {
  if (length(x) == 0L) return(numeric())
  qnorm((rank(x, ties.method = "average") - 0.5) / length(x))
}

preprocess_log <- function(message_text) {
  if (exists("pipeline_log", mode = "function")) {
    pipeline_log(message_text)
  } else {
    message(message_text)
  }
}

read_keep_samples <- function(keep_samples) {
  if (is.null(keep_samples)) return(NULL)
  if (length(keep_samples) == 1L && file.exists(keep_samples)) {
    keep_samples <- readLines(keep_samples, warn = FALSE)
  }
  keep_samples <- normalize_sample_ids(keep_samples)
  unique(keep_samples[nzchar(keep_samples)])
}

residualize_matrix <- function(M, covariate_model) {
  qr.resid(qr(covariate_model), M)
}

finite_by_row <- function(x) {
  apply(is.finite(x), 1L, all)
}

hash_numeric_column <- function(values) {
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(as.numeric(values), path, compress = FALSE)
  unname(tools::md5sum(path)[[1L]])
}

identical_numeric_columns <- function(left, right) {
  identical(as.numeric(left), as.numeric(right))
}

collapse_aligned_covariates <- function(aligned_covariates) {
  if (!length(aligned_covariates) || is.null(names(aligned_covariates))) {
    stop("Aligned covariates must be a named non-empty list.", call. = FALSE)
  }
  sample_ids <- rownames(aligned_covariates[[1L]])
  if (is.null(sample_ids) || anyDuplicated(sample_ids)) {
    stop("Aligned covariates must contain unique sample IDs.", call. = FALSE)
  }
  for (matrix in aligned_covariates) {
    if (!identical(rownames(matrix), sample_ids)) {
      stop("Aligned covariate matrices must use the same sample order.", call. = FALSE)
    }
    if (is.null(colnames(matrix)) || anyDuplicated(colnames(matrix))) {
      stop("Covariate matrices must contain unique covariate IDs.", call. = FALSE)
    }
  }

  source_rows <- do.call(rbind, lapply(names(aligned_covariates), function(modality) {
    matrix <- aligned_covariates[[modality]]
    data.frame(
      source_modality = modality,
      original_name = colnames(matrix),
      source_index = seq_len(ncol(matrix)),
      rank_before = qr(matrix)$rank,
      stringsAsFactors = FALSE
    )
  }))
  source_rows$final_name <- NA_character_
  source_rows$deduplicated_to <- NA_character_
  source_rows$md5 <- NA_character_

  source_values <- lapply(seq_len(nrow(source_rows)), function(index) {
    row <- source_rows[index, ]
    aligned_covariates[[row$source_modality]][, row$source_index]
  })
  conflicting_names <- unique(source_rows$original_name)[vapply(
    unique(source_rows$original_name),
    function(original_name) {
      indices <- which(source_rows$original_name == original_name)
      values <- source_values[indices]
      any(!vapply(values[-1L], function(candidate) {
        identical_numeric_columns(values[[1L]], candidate)
      }, logical(1L)))
    },
    logical(1L)
  )]

  output <- matrix(numeric(), nrow = length(sample_ids), ncol = 0L)
  rownames(output) <- sample_ids
  representative_indices <- integer()
  representative_names <- character()
  for (source_index in seq_len(nrow(source_rows))) {
    source_row <- source_rows[source_index, ]
    values <- source_values[[source_index]]
    duplicate_index <- which(vapply(representative_indices, function(index) {
      identical(source_rows$original_name[[index]], source_row$original_name) &&
        identical_numeric_columns(source_values[[index]], values)
    }, logical(1L)))

    if (length(duplicate_index)) {
      final_name <- representative_names[[duplicate_index[[1L]]]]
    } else {
      final_name <- if (source_row$original_name %in% conflicting_names) {
        paste(
          source_row$source_modality,
          source_row$original_name,
          sep = "::"
        )
      } else {
        source_row$original_name
      }
      output <- cbind(output, values)
      colnames(output)[[ncol(output)]] <- final_name
      representative_indices <- c(representative_indices, source_index)
      representative_names <- c(representative_names, final_name)
    }
    source_rows$final_name[[source_index]] <- final_name
    source_rows$deduplicated_to[[source_index]] <- final_name
    source_rows$md5[[source_index]] <- hash_numeric_column(values)
  }

  source_rows$rank_after <- qr(output)$rank
  list(
    matrix = output,
    provenance = data.table::as.data.table(source_rows)[, source_index := NULL]
  )
}

make_genotype_covariates <- function(modality_covariates, samples) {
  expected <- required_joint_modalities()
  if (!identical(sort(names(modality_covariates)), sort(expected))) {
    stop(
      "Covariates must contain exactly: ",
      paste(expected, collapse = ", "),
      call. = FALSE
    )
  }
  aligned <- lapply(modality_covariates[expected], function(matrix) {
    missing_samples <- setdiff(samples, rownames(matrix))
    if (length(missing_samples)) {
      stop(
        "Covariate matrix is missing aligned samples: ",
        paste(missing_samples, collapse = ", "),
        call. = FALSE
      )
    }
    matrix[samples, , drop = FALSE]
  })
  collapse_aligned_covariates(aligned)
}

input_file_checksums <- function(paths) {
  paths <- unlist(paths, use.names = TRUE)
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("Input checksum paths must have names.", call. = FALSE)
  }
  if (any(!file.exists(paths))) {
    stop("Cannot checksum an input file that does not exist.", call. = FALSE)
  }
  data.table::data.table(
    input_name = names(paths),
    path = normalizePath(paths, mustWork = TRUE),
    md5 = unname(tools::md5sum(paths))
  )
}

write_covariate_provenance <- function(provenance, path) {
  required <- c(
    "source_modality", "original_name", "final_name",
    "deduplicated_to", "rank_before", "rank_after", "md5"
  )
  if (!all(required %in% names(provenance))) {
    stop("Covariate provenance is missing required columns.", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(provenance, path, sep = "\t", quote = FALSE)
  if (!file.exists(path) || file.info(path)$size == 0) {
    stop("Failed to write covariate provenance: ", path, call. = FALSE)
  }
  invisible(path)
}

validate_joint_preprocess_inputs <- function(phenotype_data, covariates_by_modality) {
  if (is.null(phenotype_data$Y) || is.null(rownames(phenotype_data$Y))) {
    stop("The joint phenotype matrix and sample IDs are required.", call. = FALSE)
  }
  if (anyDuplicated(rownames(phenotype_data$Y))) {
    stop("The joint phenotype matrix contains duplicate sample IDs.", call. = FALSE)
  }
  if (anyDuplicated(colnames(phenotype_data$Y))) {
    stop("Selected joint outcome keys must be unique.", call. = FALSE)
  }
  modalities <- as.character(phenotype_data$modalities)
  if (length(modalities) != ncol(phenotype_data$Y)) {
    stop("Phenotype modalities must match the phenotype matrix columns.", call. = FALSE)
  }
  if (!length(modalities)) {
    stop("Prepared phenotypes must contain at least one outcome.", call. = FALSE)
  }
  unsupported <- setdiff(unique(modalities), required_joint_modalities())
  if (length(unsupported)) {
    stop(
      "Prepared phenotypes contain an unsupported modality: ",
      paste(unsupported, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(phenotype_data$metadata) != ncol(phenotype_data$Y)) {
    stop("Phenotype metadata must match the phenotype matrix columns.", call. = FALSE)
  }
  if (!identical(
    sort(names(covariates_by_modality)),
    sort(required_joint_modalities())
  )) {
    stop(
      "Covariates must contain every required joint modality.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

prepare_joint_window_data <- function(
    window,
    phenotype_data,
    dosage,
    covariates_by_modality,
    keep_samples = NULL,
    min_genotype_variance = 1e-8,
    min_phenotype_variance = 1e-8
) {
  required_window <- c("window_id", "chrom", "start", "end")
  if (!all(required_window %in% names(window))) {
    stop("Window row is missing required columns.", call. = FALSE)
  }
  if (is.null(dosage$X) || is.null(rownames(dosage$X))) {
    stop("The dosage matrix and sample IDs are required.", call. = FALSE)
  }
  if (anyDuplicated(rownames(dosage$X)) || anyDuplicated(colnames(dosage$X))) {
    stop("Dosage sample and variant IDs must be unique.", call. = FALSE)
  }
  validate_joint_preprocess_inputs(phenotype_data, covariates_by_modality)

  preprocess_log(sprintf(
    "Aligning joint samples for window %s.",
    as.character(window$window_id)
  ))
  covariate_samples <- Reduce(intersect, lapply(covariates_by_modality, rownames))
  samples <- dosage$sample_ids[
    dosage$sample_ids %in% rownames(phenotype_data$Y) &
      dosage$sample_ids %in% covariate_samples
  ]
  requested_samples <- read_keep_samples(keep_samples)
  if (!is.null(requested_samples)) {
    samples <- samples[samples %in% requested_samples]
  }
  if (!length(samples)) {
    stop("No shared samples for window: ", window$window_id, call. = FALSE)
  }
  preprocess_log(sprintf(
    "Joint sample intersection retained %d of %d genotype samples.",
    length(samples), length(dosage$sample_ids)
  ))

  X_raw <- dosage$X[samples, , drop = FALSE]
  Y_raw <- phenotype_data$Y[samples, , drop = FALSE]
  aligned_covariates <- lapply(covariates_by_modality, function(matrix) {
    matrix[samples, , drop = FALSE]
  })

  complete <- finite_by_row(X_raw) & finite_by_row(Y_raw)
  for (matrix in aligned_covariates) {
    complete <- complete & finite_by_row(matrix)
  }
  if (!all(complete)) {
    preprocess_log(sprintf(
      "Removing %d incomplete joint samples.",
      sum(!complete)
    ))
    samples <- samples[complete]
    X_raw <- X_raw[complete, , drop = FALSE]
    Y_raw <- Y_raw[complete, , drop = FALSE]
    aligned_covariates <- lapply(aligned_covariates, function(matrix) {
      matrix[complete, , drop = FALSE]
    })
  }
  if (!length(samples)) {
    stop("No complete samples for window: ", window$window_id, call. = FALSE)
  }

  phenotype_modalities <- as.character(phenotype_data$modalities)
  phenotype_metadata <- data.table::copy(phenotype_data$metadata)
  raw_phenotype_sd <- apply(Y_raw, 2L, sd)
  keep_phenotype <- is.finite(raw_phenotype_sd) &
    raw_phenotype_sd > min_phenotype_variance
  if (!any(keep_phenotype)) {
    stop("No usable phenotypes for window: ", window$window_id, call. = FALSE)
  }
  if (any(!keep_phenotype)) {
    preprocess_log(sprintf(
      "Removing %d raw zero-variance outcomes.",
      sum(!keep_phenotype)
    ))
  }
  Y_raw <- Y_raw[, keep_phenotype, drop = FALSE]
  phenotype_modalities <- phenotype_modalities[keep_phenotype]
  phenotype_metadata <- phenotype_metadata[keep_phenotype]

  preprocess_log("Applying rank-based inverse-normal transformation to each outcome.")
  Y_int <- vapply(
    seq_len(ncol(Y_raw)),
    function(index) rank_int(Y_raw[, index]),
    numeric(nrow(Y_raw))
  )
  Y_int <- as.matrix(Y_int)
  rownames(Y_int) <- samples
  colnames(Y_int) <- colnames(Y_raw)

  preprocess_log("Constructing the aligned union of joint genotype covariates.")
  genotype_union <- collapse_aligned_covariates(aligned_covariates)
  genotype_model <- cbind(genotype_union$matrix, intercept = 1)
  X_resid <- residualize_matrix(X_raw, genotype_model)

  Y_resid <- matrix(NA_real_, nrow = nrow(Y_int), ncol = ncol(Y_int))
  colnames(Y_resid) <- colnames(Y_int)
  rownames(Y_resid) <- samples
  phenotype_covariate_rank <- stats::setNames(
    integer(length(required_joint_modalities())),
    required_joint_modalities()
  )
  for (modality in required_joint_modalities()) {
    indices <- which(phenotype_modalities == modality)
    if (!length(indices)) {
      preprocess_log(sprintf(
        paste(
          "No %s outcomes were selected for this window;",
          "skipping phenotype residualization."
        ),
        modality
      ))
      next
    }
    preprocess_log(sprintf(
      "Residualizing %d %s outcomes against %d covariates.",
      length(indices), modality, ncol(aligned_covariates[[modality]])
    ))
    phenotype_model <- cbind(
      aligned_covariates[[modality]],
      intercept = 1
    )
    phenotype_qr <- qr(phenotype_model)
    Y_resid[, indices] <- qr.resid(
      phenotype_qr,
      Y_int[, indices, drop = FALSE]
    )
    phenotype_covariate_rank[[modality]] <- phenotype_qr$rank
  }

  phenotype_resid_sd <- apply(Y_resid, 2L, sd)
  keep_phenotype_resid <- is.finite(phenotype_resid_sd) &
    phenotype_resid_sd > min_phenotype_variance
  if (!any(keep_phenotype_resid)) {
    stop(
      "No usable residualized phenotypes for window: ",
      window$window_id,
      call. = FALSE
    )
  }
  if (any(!keep_phenotype_resid)) {
    preprocess_log(sprintf(
      "Removing %d zero-variance residualized outcomes.",
      sum(!keep_phenotype_resid)
    ))
  }
  Y_resid <- Y_resid[, keep_phenotype_resid, drop = FALSE]
  phenotype_modalities <- phenotype_modalities[keep_phenotype_resid]
  phenotype_metadata <- phenotype_metadata[keep_phenotype_resid]
  preprocess_log("Centering and scaling each residualized outcome to unit variance.")
  Y_scaled <- scale(Y_resid, center = TRUE, scale = TRUE)
  colnames(Y_scaled) <- colnames(Y_resid)
  rownames(Y_scaled) <- samples

  genotype_sd <- apply(X_resid, 2L, sd)
  keep_variant <- is.finite(genotype_sd) & genotype_sd > min_genotype_variance
  if (!any(keep_variant)) {
    stop("No usable variants for window: ", window$window_id, call. = FALSE)
  }
  if (any(!keep_variant)) {
    preprocess_log(sprintf(
      "Removing %d non-finite or zero-variance residualized variants.",
      sum(!keep_variant)
    ))
  }
  X_resid <- X_resid[, keep_variant, drop = FALSE]

  variant_metadata <- data.table::copy(dosage$metadata)[keep_variant]
  variant_metadata[, variant_id := dosage$variant_ids[keep_variant]]
  phenotype_metadata[, outcome_key := colnames(Y_scaled)]
  phenotype_metadata[, modality := phenotype_modalities]

  preprocess_log(sprintf(
    "Joint preprocessing complete: %d samples, %d variants, and %d outcomes.",
    length(samples), ncol(X_resid), ncol(Y_scaled)
  ))
  list(
    window = as.list(window),
    X = X_resid,
    Y = Y_scaled,
    variant_metadata = variant_metadata,
    phenotype_metadata = phenotype_metadata,
    samples = samples,
    genotype_covariates = genotype_union$matrix,
    covariate_provenance = genotype_union$provenance,
    covariate_rank = qr(genotype_model)$rank,
    phenotype_covariate_rank = phenotype_covariate_rank,
    qc = list(
      window_id = as.character(window$window_id),
      input_samples = length(dosage$sample_ids),
      shared_samples = length(samples),
      input_variants = ncol(dosage$X),
      retained_variants = ncol(X_resid),
      excluded_variants = sum(!keep_variant),
      input_phenotypes = ncol(phenotype_data$Y),
      retained_phenotypes = ncol(Y_scaled),
      excluded_phenotypes = ncol(phenotype_data$Y) - ncol(Y_scaled),
      excluded_samples = length(dosage$sample_ids) - length(samples),
      covariate_rank = qr(genotype_model)$rank,
      phenotype_covariate_rank = phenotype_covariate_rank,
      min_genotype_variance = min_genotype_variance,
      min_phenotype_variance = min_phenotype_variance,
      transformation = "rank_int_then_modality_residualization_then_unit_variance"
    )
  )
}
