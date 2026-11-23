#!/usr/bin/env Rscript

file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
scripts <- file.path(dirname(normalizePath(file)), "..", "scripts")
source(file.path(scripts, "susie_leads.R"))
source(file.path(scripts, "trans_window_prior.R"))

tmp <- tempfile("susie-test-")
dir.create(tmp)
paths <- stats::setNames(file.path(tmp, c("expression.tsv", "splicing.tsv", "protein.tsv")),
                         c("expression", "splicing", "protein"))
rows <- data.table::data.table(
  molecular_trait_id = c("expr", "expr", "expr", "expr", "other"),
  cs_id = c("L1", "L1", "L2", "L3", "L1"),
  pip = c(0.1, 0.9, 0.8, 0.7, 0.95),
  chromosome = c("chr1", "chr1", "1", "1", "1"),
  position = c(100, 101, 101, 102, 103),
  ref = "A",
  alt = "G"
)
data.table::fwrite(rows, paths[["expression"]], sep = "\t")
data.table::fwrite(rows[0], paths[["splicing"]], sep = "\t")
data.table::fwrite(rows[0], paths[["protein"]], sep = "\t")
leads <- read_susie_leads(paths)
stopifnot(nrow(leads) == 4L, leads[cs_id == "L1" & molecular_trait_id == "expr", position] == "101")

prepared <- list(
  X = matrix(0, nrow = 2L, ncol = 2L, dimnames = list(NULL, c("v1", "v2"))),
  phenotype_metadata = data.frame(modality = "expression", phenotype_id = "expr"),
  variant_metadata = data.frame(CHROM = c("1", "1"), POS = c(101, 102), REF = "A", ALT = "G")
)
matched <- match_susie_leads(leads, prepared)
stopifnot(identical(matched$variant_ids, c("v1", "v2")),
          matched$n_leads == 3L, matched$n_unmatched == 0L,
          matched$n_matched_variants == 2L)
prepared$variant_metadata$ALT[[2L]] <- "T"
matched <- match_susie_leads(leads, prepared)
stopifnot(identical(matched$variant_ids, "v1"), matched$n_unmatched == 1L)

mash_data <- list(Bhat = matrix(0, nrow = 3L, ncol = 2L,
                                dimnames = list(c("v1", "v2", "v3"), NULL)))
select_mashr_covariance_rows <- function(mash_data, n_pca, lfsr_threshold) {
  list(rows = c(2L, 3L), significant_n = 2L, fallback_used = FALSE)
}
selected <- select_covariance_training_rows(mash_data, 2L, 2L, 0.05, c("v1", "v2", "v1"))
stopifnot(identical(selected$rows, c(1L, 2L)), !selected$susie_fallback,
          selected$susie_matched_variants == 2L)
fallback <- select_covariance_training_rows(mash_data, 2L, 2L, 0.05, "v1")
stopifnot(identical(fallback$rows, c(2L, 3L)), fallback$susie_fallback,
          fallback$susie_matched_variants == 1L)
empty <- select_covariance_training_rows(mash_data, 2L, 2L, 0.05, character())
stopifnot(empty$susie_fallback, empty$susie_matched_variants == 0L)
baseline <- select_covariance_training_rows(mash_data, 2L, 2L, 0.05, NULL)
stopifnot(!baseline$susie_mode, !baseline$susie_fallback,
          identical(baseline$rows, c(2L, 3L)))
stopifnot(inherits(try(select_covariance_training_rows(mash_data, 2L, 2L, 0.05, "missing"),
                       silent = TRUE), "try-error"))
unlink(paths)
unlink(tmp)
cat("SuSiE lead selection tests passed.\n")
