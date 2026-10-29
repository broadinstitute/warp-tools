module_file <- normalizePath(sys.frame(1L)$ofile)
module_dir <- dirname(module_file)
source(file.path(module_dir, "trans_window_logging.R"))
source(file.path(module_dir, "trans_window_prior.R"))

validate_positive_integer <- function(value, label) {
  if (
    length(value) != 1L || !is.numeric(value) || is.na(value) ||
    !is.finite(value) || value < 1L || value != as.integer(value)
  ) {
    stop(label, " must be a positive integer.", call. = FALSE)
  }
  as.integer(value)
}

make_model_config <- function(
    start_L = 10L,
    step_L = 5L,
    max_L = 40L,
    greedy_lbf_cutoff = 1,
    max_iter = 100L,
    tol = 1e-4,
    coverage = 0.95,
    min_abs_corr = 0.5,
    n_thread = 1L,
    mashr_n_pca = 5L,
    mashr_seed = NULL,
    mashr_strong_lfsr = 0.05
) {
  start_L <- validate_positive_integer(start_L, "start_L")
  step_L <- validate_positive_integer(step_L, "step_L")
  max_L <- validate_positive_integer(max_L, "max_L")
  if (start_L > max_L) {
    stop("start_L cannot exceed max_L.", call. = FALSE)
  }
  if ((max_L - start_L) %% step_L != 0L) {
    stop("The greedy schedule must reach max_L exactly.", call. = FALSE)
  }
  if (
    length(greedy_lbf_cutoff) != 1L || is.na(greedy_lbf_cutoff) ||
    !is.finite(greedy_lbf_cutoff)
  ) {
    stop("greedy_lbf_cutoff must be one finite number.", call. = FALSE)
  }
  max_iter <- validate_positive_integer(max_iter, "max_iter")
  n_thread <- validate_positive_integer(n_thread, "n_thread")
  mashr_n_pca <- validate_positive_integer(mashr_n_pca, "mashr_n_pca")
  if (
    length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0 ||
    length(coverage) != 1L || is.na(coverage) ||
    coverage <= 0 || coverage >= 1 ||
    length(min_abs_corr) != 1L || is.na(min_abs_corr) ||
    min_abs_corr < 0 || min_abs_corr > 1 ||
    length(mashr_strong_lfsr) != 1L || is.na(mashr_strong_lfsr) ||
    mashr_strong_lfsr <= 0 || mashr_strong_lfsr >= 1
  ) {
    stop("Model numeric thresholds are outside their valid ranges.", call. = FALSE)
  }
  if (!is.null(mashr_seed)) {
    mashr_seed <- as.integer(mashr_seed)
    if (is.na(mashr_seed)) stop("mashr_seed must be an integer.", call. = FALSE)
  }
  list(
    start_L = start_L,
    step_L = step_L,
    max_L = max_L,
    greedy_lbf_cutoff = as.numeric(greedy_lbf_cutoff),
    max_iter = max_iter,
    tol = as.numeric(tol),
    coverage = as.numeric(coverage),
    min_abs_corr = as.numeric(min_abs_corr),
    n_thread = n_thread,
    mashr_n_pca = mashr_n_pca,
    mashr_seed = mashr_seed,
    mashr_strong_lfsr = as.numeric(mashr_strong_lfsr)
  )
}

validate_prior_for_outcomes <- function(prior, n_outcomes) {
  if (!inherits(prior, "mash_prior") || !length(prior$xUlist)) {
    stop("The mvSuSiE prior must be a non-empty mash prior.", call. = FALSE)
  }
  validate_mashr_covariances(prior$xUlist, n_outcomes)
  if (
    length(prior$pi) != length(prior$xUlist) ||
    any(!is.finite(prior$pi)) || any(prior$pi < 0) ||
    abs(sum(prior$pi) - 1) > 1e-8
  ) {
    stop("The mvSuSiE prior has invalid fixed mixture weights.", call. = FALSE)
  }
  invisible(TRUE)
}

supported_component_count <- function(single_effect_lfsr, fitted_L) {
  if (is.null(single_effect_lfsr) || all(is.na(single_effect_lfsr))) return(0L)
  if (is.matrix(single_effect_lfsr) && nrow(single_effect_lfsr) == fitted_L) {
    return(sum(apply(single_effect_lfsr < 0.05, 1L, any, na.rm = TRUE)))
  }
  if (length(single_effect_lfsr) == fitted_L) {
    return(sum(single_effect_lfsr < 0.05, na.rm = TRUE))
  }
  0L
}

validate_mvsusie_round <- function(fit, requested_L) {
  if (!isTRUE(fit$converged)) {
    stop("The mvSuSiE greedy round did not converge at L = ", requested_L, ".", call. = FALSE)
  }
  if (
    is.null(fit$alpha) || is.null(fit$lbf) ||
    any(!is.finite(fit$alpha)) || any(!is.finite(fit$lbf))
  ) {
    stop("The mvSuSiE greedy round returned non-finite values.", call. = FALSE)
  }
  invisible(TRUE)
}

fit_mvsusie_greedy_schedule <- function(
    X,
    Y,
    prior,
    start_L = 10L,
    step_L = 5L,
    max_L = 40L,
    greedy_lbf_cutoff = 1,
    fit_fun = mvsusieR::mvsusie,
    ...
) {
  config <- make_model_config(
    start_L = start_L,
    step_L = step_L,
    max_L = max_L,
    greedy_lbf_cutoff = greedy_lbf_cutoff
  )
  validate_prior_for_outcomes(prior, ncol(Y))

  requested_values <- seq.int(
    config$start_L,
    config$max_L,
    by = config$step_L
  )
  previous_fit <- NULL
  history_rows <- vector("list", length(requested_values))
  completed_rounds <- 0L

  for (round_index in seq_along(requested_values)) {
    requested_L <- requested_values[[round_index]]
    pipeline_log(sprintf(
      "Starting greedy mvSuSiE round %d at L = %d.",
      round_index, requested_L
    ))
    fit <- fit_fun(
      X = X,
      Y = Y,
      L = requested_L,
      prior_variance = prior,
      model_init = previous_fit,
      ...
    )
    validate_mvsusie_round(fit, requested_L)

    fitted_L <- as.integer(nrow(fit$alpha))
    minimum_lbf <- min(as.numeric(fit$lbf))
    saturated <- minimum_lbf < config$greedy_lbf_cutoff
    at_maximum <- requested_L == config$max_L
    action <- if (saturated) {
      "saturated"
    } else if (at_maximum) {
      "maximum"
    } else {
      "continue"
    }
    credible_set_count <- if (is.null(fit$sets$cs)) 0L else length(fit$sets$cs)
    history_rows[[round_index]] <- data.table::data.table(
      round = as.integer(round_index),
      requested_L = as.integer(requested_L),
      fitted_L = fitted_L,
      niter = as.integer(fit$niter),
      minimum_lbf = minimum_lbf,
      credible_set_count = as.integer(credible_set_count),
      supported_component_count = supported_component_count(
        fit$single_effect_lfsr,
        fitted_L
      ),
      maximum_alpha = max(as.numeric(fit$alpha)),
      action = action
    )
    completed_rounds <- round_index
    pipeline_log(sprintf(
      paste(
        "Greedy round %d complete: fitted L=%d, min lbf=%.6g,",
        "action=%s."
      ),
      round_index, fitted_L, minimum_lbf, action
    ))
    previous_fit <- fit
    if (saturated || at_maximum) break
  }

  list(
    fit = previous_fit,
    history = data.table::rbindlist(history_rows[seq_len(completed_rounds)])
  )
}

fit_window_mvsusie <- function(prepared, config) {
  if (
    !is.matrix(prepared$X) || !is.matrix(prepared$Y) ||
    nrow(prepared$X) != nrow(prepared$Y)
  ) {
    stop("Prepared genotype and phenotype matrices have incompatible dimensions.", call. = FALSE)
  }
  if (any(!is.finite(prepared$X)) || any(!is.finite(prepared$Y))) {
    stop("Prepared genotype and phenotype matrices must be finite.", call. = FALSE)
  }

  pipeline_log("Computing all-SNP associations for the joint mashr prior.")
  marginal <- compute_marginal_bhat_shat_matrix(prepared$X, prepared$Y)
  mashr_training <- learn_joint_mashr_prior(
    Bhat = marginal$Bhat,
    Shat = marginal$Shat,
    n_pca = config$mashr_n_pca,
    seed = config$mashr_seed,
    strong_lfsr = config$mashr_strong_lfsr
  )
  raw_prior <- mashr_training$raw_prior
  validate_prior_for_outcomes(raw_prior, ncol(prepared$Y))
  raw_covariance_values <- unlist(raw_prior$xUlist, use.names = FALSE)
  prepared_prior <- prepare_mashr_prior_for_mvsusie(raw_prior, prepared$Y)
  validate_prior_for_outcomes(prepared_prior, ncol(prepared$Y))
  converted_covariance_values <- unlist(
    prepared_prior$xUlist,
    use.names = FALSE
  )
  outcome_se <- attr(prepared_prior, "mvsusie_outcome_se_scale")
  mashr_training$converted_prior <- prepared_prior
  mashr_training$outcome_se <- outcome_se
  mashr_training$raw_covariance_range <- range(raw_covariance_values)
  mashr_training$converted_covariance_range <- range(converted_covariance_values)

  residual_variance <- stats::cov(prepared$Y)
  if (any(!is.finite(residual_variance))) {
    stop("The initial residual covariance contains non-finite values.", call. = FALSE)
  }
  pipeline_log("Using fixed mashr weights and the fixed initial residual covariance.")
  pipeline_log("Starting mvSuSiE with verbose iteration output.")
  scheduled <- fit_mvsusie_greedy_schedule(
    X = prepared$X,
    Y = prepared$Y,
    prior = prepared_prior,
    start_L = config$start_L,
    step_L = config$step_L,
    max_L = config$max_L,
    greedy_lbf_cutoff = config$greedy_lbf_cutoff,
    residual_variance = residual_variance,
    standardize = TRUE,
    intercept = FALSE,
    estimate_residual_variance = FALSE,
    estimate_prior_variance = FALSE,
    estimate_prior_mixture_weights = FALSE,
    coverage = config$coverage,
    min_abs_corr = config$min_abs_corr,
    precompute_cache = TRUE,
    n_thread = config$n_thread,
    max_iter = config$max_iter,
    tol = config$tol,
    verbose = TRUE
  )
  fit <- scheduled$fit
  pipeline_log(sprintf(
    "Final mvSuSiE fit completed at L = %d after %d iterations.",
    nrow(fit$alpha), fit$niter
  ))

  list(
    fit = fit,
    greedy_history = scheduled$history,
    mashr_training = mashr_training,
    metadata = list(
      window_id = prepared$qc$window_id,
      prior = "mashr_pca_only",
      prior_components = length(prepared_prior$xUlist),
      pca_requested = mashr_training$pca_requested,
      pca_used = mashr_training$pca_used,
      pca_returned = mashr_training$pca_returned,
      mash_model_training_scope = mashr_training$mash_model_training_scope,
      mash_model_training_n = mashr_training$mash_model_training_n,
      covariance_training_scope = mashr_training$covariance_training_scope,
      covariance_training_n = mashr_training$covariance_training_n,
      covariance_significant_n = mashr_training$covariance_significant_n,
      covariance_selection_lfsr = mashr_training$covariance_selection_lfsr,
      covariance_selection_fallback_used =
        mashr_training$covariance_selection_fallback_used,
      covariance_input_method = mashr_training$covariance_input_method,
      prior_mixture_weights_mode = "fixed_from_mashr",
      prior_scale_conversion = "preserve_mashr_effect_covariance",
      prior_outcome_se_min = min(outcome_se),
      prior_outcome_se_max = max(outcome_se),
      residual_variance_mode = "fixed_initial_covariance",
      mvsusieR_version = as.character(utils::packageVersion("mvsusieR")),
      config = config,
      L_final = as.integer(nrow(fit$alpha)),
      converged = isTRUE(fit$converged),
      niter = fit$niter
    )
  )
}

empty_credible_set_summary <- function() {
  data.table::data.table(
    component = integer(),
    credible_set_size = integer(),
    sentinel_variant_id = character(),
    sentinel_alpha = numeric(),
    coverage = numeric(),
    purity_min = numeric(),
    purity_mean = numeric()
  )
}

empty_credible_set_members <- function() {
  data.table::data.table(
    component = integer(),
    variant_id = character(),
    alpha = numeric(),
    pip = numeric(),
    is_sentinel = logical()
  )
}

extract_variant_pips <- function(fit, prepared) {
  data.table::data.table(
    variant_id = colnames(prepared$X),
    pip = as.numeric(fit$pip)
  )
}

finite_summary <- function(values, fun) {
  values <- values[is.finite(values)]
  if (!length(values)) return(NA_real_)
  fun(values)
}

extract_credible_set_tables <- function(fit, prepared, config) {
  cs_obj <- susieR::susie_get_cs(
    fit,
    X = prepared$X,
    coverage = config$coverage,
    min_abs_corr = config$min_abs_corr
  )
  if (!length(cs_obj$cs)) {
    return(list(
      summary = empty_credible_set_summary(),
      members = empty_credible_set_members()
    ))
  }
  purity <- cs_obj$purity
  summary_rows <- vector("list", length(cs_obj$cs))
  member_rows <- vector("list", length(cs_obj$cs))
  for (i in seq_along(cs_obj$cs)) {
    members <- cs_obj$cs[[i]]
    component <- as.integer(sub("^L", "", names(cs_obj$cs)[[i]]))
    purity_row <- if (is.null(dim(purity))) purity else purity[i, ]
    purity_values <- suppressWarnings(as.numeric(unlist(
      purity_row,
      use.names = FALSE
    )))
    alpha <- as.numeric(fit$alpha[component, members])
    sentinel_index <- which.max(alpha)
    member_rows[[i]] <- data.table::data.table(
      component = component,
      variant_id = colnames(prepared$X)[members],
      alpha = alpha,
      pip = as.numeric(fit$pip[members]),
      is_sentinel = seq_along(members) == sentinel_index
    )
    summary_rows[[i]] <- data.table::data.table(
      component = component,
      credible_set_size = as.integer(length(members)),
      sentinel_variant_id = colnames(prepared$X)[members[[sentinel_index]]],
      sentinel_alpha = alpha[[sentinel_index]],
      coverage = config$coverage,
      purity_min = finite_summary(purity_values, min),
      purity_mean = finite_summary(purity_values, mean)
    )
  }
  list(
    summary = data.table::rbindlist(summary_rows, fill = TRUE),
    members = data.table::rbindlist(member_rows, fill = TRUE)
  )
}

get_fit_mu2 <- function(fit) {
  if (!is.null(fit$mu2_diag)) return(fit$mu2_diag)
  if (!is.null(fit$mu2)) return(fit$mu2)
  stop("The mvSuSiE fit contains neither mu2_diag nor mu2.", call. = FALSE)
}

as_component_outcome_matrix <- function(value, n_component, n_outcome, name) {
  if (is.null(value) || (length(value) == 1L && is.na(value))) {
    return(matrix(NA_real_, nrow = n_component, ncol = n_outcome))
  }
  value <- as.matrix(value)
  if (!identical(dim(value), c(n_component, n_outcome))) {
    stop(name, " must have one row per component and one column per outcome.", call. = FALSE)
  }
  value
}

extract_component_feature_support <- function(fit, prepared) {
  n_component <- nrow(fit$alpha)
  n_outcome <- ncol(prepared$Y)
  lfsr <- as_component_outcome_matrix(
    fit$single_effect_lfsr, n_component, n_outcome, "single_effect_lfsr"
  )
  outcome_lbf <- as_component_outcome_matrix(
    fit$lbf_outcome, n_component, n_outcome, "lbf_outcome"
  )
  idx <- expand.grid(
    component = seq_len(n_component),
    outcome_index = seq_len(n_outcome)
  )
  metadata <- data.table::as.data.table(prepared$phenotype_metadata)
  metadata <- metadata[match(colnames(prepared$Y), outcome_key)]
  if (anyNA(metadata$outcome_key)) {
    stop("Phenotype metadata does not match the prepared outcome matrix.", call. = FALSE)
  }
  data.table::data.table(
    component = idx$component,
    outcome_key = metadata$outcome_key[idx$outcome_index],
    modality = metadata$modality[idx$outcome_index],
    phenotype_id = metadata$phenotype_id[idx$outcome_index],
    single_effect_lfsr = as.vector(lfsr),
    outcome_lbf = as.vector(outcome_lbf)
  )
}
