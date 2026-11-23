module_file <- normalizePath(sys.frame(1L)$ofile)
module_dir <- dirname(module_file)
source(file.path(module_dir, "trans_window_logging.R"))

validate_marginal_summary_statistics <- function(Bhat, Shat) {
  if (!is.matrix(Bhat) || !is.numeric(Bhat)) {
    stop("Bhat must be a numeric matrix.", call. = FALSE)
  }
  if (!is.matrix(Shat) || !is.numeric(Shat)) {
    stop("Shat must be a numeric matrix.", call. = FALSE)
  }
  if (!identical(dim(Bhat), dim(Shat))) {
    stop("Bhat and Shat must have identical dimensions.", call. = FALSE)
  }
  if (any(!is.finite(Bhat)) || any(!is.finite(Shat))) {
    stop("Bhat and Shat must contain only finite values.", call. = FALSE)
  }
  if (any(Shat <= 0)) {
    stop("Shat must contain strictly positive values.", call. = FALSE)
  }
  invisible(TRUE)
}

make_mashr_data <- function(Bhat, Shat) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("The mashr package is required to learn a data-driven prior.", call. = FALSE)
  }
  mashr::mash_set_data(Bhat = Bhat, Shat = Shat, alpha = 0)
}

prepare_mashr_prior_for_mvsusie <- function(prior, Y) {
  if (!inherits(prior, "mash_prior") || !length(prior$xUlist)) {
    stop("prior must be a non-empty mash prior.", call. = FALSE)
  }
  if (!is.matrix(Y) || !is.numeric(Y) || ncol(Y) < 1L) {
    stop("Y must be a numeric matrix with at least one outcome.", call. = FALSE)
  }
  validate_mashr_covariances(prior$xUlist, ncol(Y))
  outcome_n <- colSums(is.finite(Y))
  outcome_sd <- apply(Y, 2L, stats::sd, na.rm = TRUE)
  outcome_se <- outcome_sd / sqrt(outcome_n)
  if (any(!is.finite(outcome_se)) || any(outcome_se <= 0)) {
    stop("Each outcome must have a finite positive standard error scale.", call. = FALSE)
  }
  automatic_scale <- tcrossprod(outcome_se)
  prior$xUlist <- lapply(prior$xUlist, function(U) U / automatic_scale)
  validate_mashr_covariances(prior$xUlist, ncol(Y))
  attr(prior, "mvsusie_outcome_se_scale") <- outcome_se
  prior
}

compute_marginal_bhat_shat_matrix <- function(X, Y, block_size = 1000L) {
  if (!is.matrix(X) || !is.numeric(X) || !is.matrix(Y) || !is.numeric(Y)) {
    stop("X and Y must be numeric matrices.", call. = FALSE)
  }
  if (nrow(X) != nrow(Y) || nrow(X) < 2L) {
    stop("X and Y must have the same number of at least two rows.", call. = FALSE)
  }
  if (any(!is.finite(X)) || any(!is.finite(Y))) {
    stop("X and Y must contain only finite values.", call. = FALSE)
  }
  if (
    length(block_size) != 1L || is.na(block_size) ||
    block_size < 1L || block_size != as.integer(block_size)
  ) {
    stop("block_size must be a positive integer.", call. = FALSE)
  }

  start_time <- proc.time()[["elapsed"]]
  n <- nrow(X)
  n_variants <- ncol(X)
  n_outcomes <- ncol(Y)
  pipeline_log(sprintf(
    "Starting marginal associations: N=%d, J=%d, R=%d.",
    n, n_variants, n_outcomes
  ))
  pipeline_log("Centering outcomes and computing their sums of squares.")
  Y_centered <- sweep(Y, 2L, colMeans(Y), "-")
  y_ss <- colSums(Y_centered^2)
  predictor_weight <- n - 1
  Bhat <- matrix(NA_real_, nrow = n_variants, ncol = n_outcomes)
  Shat <- matrix(NA_real_, nrow = n_variants, ncol = n_outcomes)
  if (!is.null(colnames(X))) rownames(Bhat) <- rownames(Shat) <- colnames(X)
  if (!is.null(colnames(Y))) colnames(Bhat) <- colnames(Shat) <- colnames(Y)
  block_starts <- seq.int(1L, n_variants, by = as.integer(block_size))

  for (block_index in seq_along(block_starts)) {
    first <- block_starts[[block_index]]
    last <- min(first + as.integer(block_size) - 1L, n_variants)
    indices <- first:last
    pipeline_log(sprintf(
      "Computing the all-SNP cross-product block %d of %d (%d-%d).",
      block_index, length(block_starts), first, last
    ))
    X_block <- X[, indices, drop = FALSE]
    X_block <- sweep(X_block, 2L, colMeans(X_block), "-")
    x_ss <- colSums(X_block^2)
    if (any(!is.finite(x_ss)) || any(x_ss <= 0)) {
      stop("X contains a non-finite or zero-variance column.", call. = FALSE)
    }
    x_scale <- sqrt(x_ss / predictor_weight)
    X_standardized <- sweep(X_block, 2L, x_scale, "/")
    block_Bhat <- crossprod(X_standardized, Y) / predictor_weight
    Bhat[indices, ] <- block_Bhat

    X_standardized_centered <- sweep(
      X_standardized,
      2L,
      colMeans(X_standardized),
      "-"
    )
    residual_crossproduct <- crossprod(X_standardized_centered, Y_centered)
    standardized_x_ss <- colSums(X_standardized_centered^2)
    residual_ss <- matrix(
      y_ss,
      nrow = length(indices),
      ncol = n_outcomes,
      byrow = TRUE
    ) - 2 * block_Bhat * residual_crossproduct +
      block_Bhat^2 * standardized_x_ss

    y_ss_matrix <- matrix(
      y_ss,
      nrow = length(indices),
      ncol = n_outcomes,
      byrow = TRUE
    )
    unstable <- !is.finite(residual_ss) |
      residual_ss <= sqrt(.Machine$double.eps) * y_ss_matrix
    if (any(unstable)) {
      pipeline_log(sprintf(
        paste(
          "Recomputing %d near-perfect association residual sums of squares",
          "directly."
        ),
        sum(unstable)
      ))
      fallback_chunk_size <- 100L
      for (outcome_index in which(colSums(unstable) > 0L)) {
        local_indices <- which(unstable[, outcome_index])
        fallback_starts <- seq.int(
          1L,
          length(local_indices),
          by = fallback_chunk_size
        )
        for (fallback_first in fallback_starts) {
          fallback_last <- min(
            fallback_first + fallback_chunk_size - 1L,
            length(local_indices)
          )
          selected <- local_indices[fallback_first:fallback_last]
          fitted_values <- sweep(
            X_standardized_centered[, selected, drop = FALSE],
            2L,
            block_Bhat[selected, outcome_index],
            "*"
          )
          residuals <- -sweep(
            fitted_values,
            1L,
            Y_centered[, outcome_index],
            "-"
          )
          residual_ss[selected, outcome_index] <- colSums(residuals^2)
        }
      }
    }
    residual_variance <- pmax(residual_ss / predictor_weight, 1e-64)
    Shat[indices, ] <- sqrt(residual_variance) / sqrt(predictor_weight)
  }
  validate_marginal_summary_statistics(Bhat, Shat)

  pipeline_log(sprintf(
    "Marginal associations complete in %.2f seconds.",
    proc.time()[["elapsed"]] - start_time
  ))
  list(Bhat = Bhat, Shat = Shat)
}

select_mashr_covariance_rows <- function(mash_data, n_pca, lfsr_threshold) {
  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf(
    "Starting the one-by-one mashr fit for strong-row selection at lfsr <= %.3g.",
    lfsr_threshold
  ))
  one_by_one <- mashr::mash_1by1(mash_data)
  strong_rows <- mashr::get_significant_results(
    one_by_one,
    thresh = lfsr_threshold
  )
  significant_n <- length(strong_rows)
  fallback_used <- FALSE
  if (significant_n < n_pca) {
    fallback_used <- TRUE
    minimum_lfsr <- apply(one_by_one$result$lfsr, 1L, min)
    strong_rows <- head(order(minimum_lfsr), n_pca)
    pipeline_log(sprintf(
      paste(
        "Only %d rows passed the strong-row threshold; using the %d rows",
        "with the smallest lfsr values for covariance training."
      ),
      significant_n, length(strong_rows)
    ))
  }
  pipeline_log(sprintf(
    paste(
      "Selected %d strong SNP rows from %d total rows in %.2f seconds",
      "for PCA covariance learning."
    ),
    length(strong_rows), nrow(mash_data$Bhat),
    proc.time()[["elapsed"]] - stage_time
  ))
  list(
    rows = as.integer(strong_rows),
    significant_n = significant_n,
    fallback_used = fallback_used
  )
}

select_covariance_training_rows <- function(
    mash_data, n_pca, pca_requested, strong_lfsr, susie_variant_ids
) {
  susie_mode <- !is.null(susie_variant_ids)
  susie_rows <- if (susie_mode) {
    match(unique(susie_variant_ids), rownames(mash_data$Bhat))
  } else {
    integer()
  }
  if (anyNA(susie_rows)) stop("SuSiE lead IDs do not match marginal statistics.", call. = FALSE)
  susie_fallback <- susie_mode && length(susie_rows) < pca_requested
  if (susie_mode && !susie_fallback) {
    selection <- list(rows = susie_rows, significant_n = NA_integer_, fallback_used = FALSE)
    pipeline_log(sprintf("Using %d unique SuSiE lead variants for covariance training.", length(susie_rows)))
  } else {
    if (susie_fallback) {
      pipeline_log(sprintf(
        "Only %d unique SuSiE leads matched; falling back to strong-association selection.",
        length(susie_rows)
      ))
    }
    selection <- select_mashr_covariance_rows(mash_data, n_pca, strong_lfsr)
  }
  selection$susie_mode <- susie_mode
  selection$susie_fallback <- susie_fallback
  selection$susie_matched_variants <- if (susie_mode) length(susie_rows) else NA_integer_
  selection
}

validate_mashr_covariances <- function(covariances, n_outcomes) {
  if (!is.list(covariances) || !length(covariances)) {
    stop("Mash covariance inputs must be a non-empty list.", call. = FALSE)
  }
  expected <- c(as.integer(n_outcomes), as.integer(n_outcomes))
  for (covariance in covariances) {
    if (!is.matrix(covariance) || !identical(dim(covariance), expected)) {
      stop("Mash covariance dimensions must match the joint outcomes.", call. = FALSE)
    }
    if (any(!is.finite(covariance))) {
      stop("Mash covariance matrices contain non-finite values.", call. = FALSE)
    }
  }
  invisible(TRUE)
}

learn_joint_mashr_prior <- function(
    Bhat,
    Shat,
    n_pca = 5L,
    seed = NULL,
    strong_lfsr = 0.05,
    cov_pca_fun = mashr::cov_pca,
    include_canonical_priors = FALSE,
    use_mashr_null_weight = FALSE,
    susie_variant_ids = NULL
) {
  validate_marginal_summary_statistics(Bhat, Shat)
  if (
    length(n_pca) != 1L || is.na(n_pca) ||
    n_pca < 1L || n_pca != as.integer(n_pca)
  ) {
    stop("n_pca must be a positive integer.", call. = FALSE)
  }
  if (
    length(strong_lfsr) != 1L || is.na(strong_lfsr) ||
    strong_lfsr <= 0 || strong_lfsr >= 1
  ) {
    stop("strong_lfsr must be between zero and one.", call. = FALSE)
  }
  if (
    length(include_canonical_priors) != 1L ||
    is.na(include_canonical_priors) ||
    !is.logical(include_canonical_priors) ||
    length(use_mashr_null_weight) != 1L ||
    is.na(use_mashr_null_weight) ||
    !is.logical(use_mashr_null_weight)
  ) {
    stop("Prior options must be single logical values.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)

  mash_data <- make_mashr_data(Bhat, Shat)
  pca_requested <- as.integer(n_pca)
  pca_used <- min(pca_requested, nrow(Bhat), ncol(Bhat))
  if (pca_used < pca_requested) {
    pipeline_log(sprintf(
      "Reducing PCA covariance inputs from %d to %d for this window.",
      pca_requested, pca_used
    ))
  }
  covariance_selection <- select_covariance_training_rows(
    mash_data, pca_used, pca_requested, strong_lfsr, susie_variant_ids
  )
  covariance_rows <- covariance_selection$rows
  stage_time <- proc.time()[["elapsed"]]
  covariance_input_method <- "pca_only"
  if (ncol(Bhat) == 1L) {
    pipeline_log(
      "Using the univariate equivalent of the PCA total covariance."
    )
    covariance_value <- mean(Bhat[covariance_rows, 1L]^2)
    covariance_value <- max(covariance_value, .Machine$double.eps)
    covariance <- matrix(
      covariance_value,
      nrow = 1L,
      ncol = 1L,
      dimnames = list(colnames(Bhat), colnames(Bhat))
    )
    pca_covariances <- list(univariate_pca_equivalent = covariance)
    covariance_input_method <- "univariate_pca_equivalent"
  } else {
    pipeline_log(sprintf(
      "Starting %d PCA covariance inputs on %d selected rows.",
      pca_used, length(covariance_rows)
    ))
    pca_covariances <- cov_pca_fun(
      mash_data,
      npc = pca_used,
      subset = covariance_rows
    )
  }
  covariance_inputs <- pca_covariances
  canonical_covariances <- list()
  if (include_canonical_priors && ncol(Bhat) > 1L) {
    pipeline_log("Adding canonical cross-outcome covariance inputs.")
    canonical_covariances <- mashr::cov_canonical(mash_data)
    names(canonical_covariances) <- paste0(
      "canonical_",
      names(canonical_covariances)
    )
    names(pca_covariances) <- paste0("pca_", names(pca_covariances))
    covariance_inputs <- c(canonical_covariances, pca_covariances)
    covariance_input_method <- "canonical_and_pca"
  } else if (include_canonical_priors) {
    pipeline_log("Skipping canonical covariances for a single outcome.")
    covariance_input_method <- "univariate_pca_equivalent"
  }
  pipeline_log(sprintf(
    "Covariance learning returned %d canonical and %d PCA matrices in %.2f seconds.",
    length(canonical_covariances),
    length(pca_covariances),
    proc.time()[["elapsed"]] - stage_time
  ))
  validate_mashr_covariances(covariance_inputs, ncol(Bhat))

  stage_time <- proc.time()[["elapsed"]]
  pipeline_log(sprintf(
    "Starting the mashr mixture fit on all %d SNP rows using %d covariance inputs.",
    nrow(Bhat), length(covariance_inputs)
  ))
  mash_fit <- mashr::mash(
    data = mash_data,
    Ulist = covariance_inputs,
    usepointmass = TRUE,
    outputlevel = 0,
    verbose = TRUE
  )
  pipeline_log(sprintf(
    "The mashr mixture fit completed in %.2f seconds.",
    proc.time()[["elapsed"]] - stage_time
  ))
  mashr_null_weight <- if (isTRUE(mash_fit$fitted_g$usepointmass)) {
    as.numeric(mash_fit$fitted_g$pi[[1L]])
  } else {
    0
  }
  requested_null_weight <- if (use_mashr_null_weight) {
    NULL
  } else {
    0
  }
  prior <- mvsusieR::create_mixture_prior(
    fitted_g = mash_fit$fitted_g,
    null_weight = requested_null_weight
  )
  fallback_to_input_covariances <- FALSE
  if (!length(prior$xUlist)) {
    prior <- mvsusieR::create_mixture_prior(
      mixture_prior = list(matrices = covariance_inputs),
      null_weight = if (use_mashr_null_weight) mashr_null_weight else 0,
      weights_tol = 0
    )
    fallback_to_input_covariances <- TRUE
  }

  validate_mashr_covariances(prior$xUlist, ncol(Bhat))
  fitted_weights <- as.numeric(prior$pi)
  if (
    length(fitted_weights) != length(prior$xUlist) ||
    any(!is.finite(fitted_weights)) || any(fitted_weights < 0) ||
    abs(sum(fitted_weights) - 1) > 1e-8
  ) {
    stop("The fitted mashr weights are invalid.", call. = FALSE)
  }

  prior_method <- if (
    !include_canonical_priors && !use_mashr_null_weight
  ) {
    "mashr_pca_only"
  } else {
    paste0(
      "mashr_", covariance_input_method, "_",
      if (use_mashr_null_weight) "estimated_null" else "null_weight_zero"
    )
  }

  list(
    raw_prior = prior,
    Bhat = Bhat,
    Shat = Shat,
    mash_model_training_scope = "all_snps_in_window",
    mash_model_training_n = nrow(Bhat),
    covariance_training_scope = if (covariance_selection$susie_mode && !covariance_selection$susie_fallback) {
      "susie_lead_snps_in_window"
    } else {
      "strong_snps_in_window"
    },
    covariance_training_n = length(covariance_rows),
    covariance_significant_n = covariance_selection$significant_n,
    covariance_selection_lfsr = if (covariance_selection$susie_mode && !covariance_selection$susie_fallback) NA_real_ else strong_lfsr,
    covariance_selection_fallback_used = covariance_selection$fallback_used,
    covariance_training_source = if (covariance_selection$susie_mode && !covariance_selection$susie_fallback) "susie_leads" else "strong_associations",
    susie_lead_mode = covariance_selection$susie_mode,
    susie_matched_variants = covariance_selection$susie_matched_variants,
    susie_lead_fallback = covariance_selection$susie_fallback,
    covariance_rows = rownames(Bhat)[covariance_rows],
    pca_requested = pca_requested,
    pca_used = pca_used,
    pca_returned = length(pca_covariances),
    canonical_returned = length(canonical_covariances),
    covariance_input_method = covariance_input_method,
    n_covariance_inputs = length(covariance_inputs),
    n_prior_components = length(prior$xUlist),
    fallback_to_input_covariances = fallback_to_input_covariances,
    prior_method = prior_method,
    use_mashr_null_weight = use_mashr_null_weight,
    mashr_estimated_null_weight = mashr_null_weight,
    mvsusie_null_weight = prior$null_weight,
    fitted_g = mash_fit$fitted_g,
    fitted_weights = fitted_weights,
    seed = seed
  )
}
