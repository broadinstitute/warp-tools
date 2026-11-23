make_mvsusie_prior_learner <- function(
    include_canonical_priors = FALSE,
    use_mashr_null_weight = FALSE,
    susie_variant_ids = NULL
) {
  if (
    length(include_canonical_priors) != 1L ||
    !is.logical(include_canonical_priors) ||
    is.na(include_canonical_priors) ||
    length(use_mashr_null_weight) != 1L ||
    !is.logical(use_mashr_null_weight) ||
    is.na(use_mashr_null_weight)
  ) {
    stop("Prior options must be single logical values.", call. = FALSE)
  }

  function(
      Bhat,
      Shat,
      n_pca = 5L,
      seed = NULL,
      strong_lfsr = 0.05,
      cov_pca_fun = mashr::cov_pca
  ) {
    learn_joint_mashr_prior(
      Bhat = Bhat,
      Shat = Shat,
      n_pca = n_pca,
      seed = seed,
      strong_lfsr = strong_lfsr,
      cov_pca_fun = cov_pca_fun,
      include_canonical_priors = include_canonical_priors,
      use_mashr_null_weight = use_mashr_null_weight,
      susie_variant_ids = susie_variant_ids
    )
  }
}
