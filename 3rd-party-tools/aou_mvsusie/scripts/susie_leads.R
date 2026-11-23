read_susie_leads <- function(paths) {
  modalities <- c("expression", "splicing", "protein")
  if (!identical(names(paths), modalities) || any(!nzchar(paths))) {
    stop("Provide one annotated SuSiE TSV for each modality.", call. = FALSE)
  }
  required <- c("molecular_trait_id", "cs_id", "pip", "chromosome", "position", "ref", "alt")
  leads <- lapply(modalities, function(modality) {
    path <- paths[[modality]]
    if (!file.exists(path)) stop("SuSiE TSV not found: ", path, call. = FALSE)
    header <- names(data.table::fread(path, nrows = 0L))
    missing <- setdiff(required, header)
    if (length(missing)) {
      stop(modality, " SuSiE TSV lacks columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    dt <- data.table::fread(
      path, select = required,
      colClasses = list(character = setdiff(required, "pip"))
    )
    dt[, pip := suppressWarnings(as.numeric(pip))]
    if (anyNA(dt[, ..required]) || any(!is.finite(dt$pip)) ||
        any(dt$pip < 0 | dt$pip > 1) ||
        any(!nzchar(trimws(dt$molecular_trait_id))) ||
        any(!nzchar(trimws(dt$cs_id))) ||
        any(!nzchar(trimws(dt$chromosome))) ||
        any(!grepl("^[0-9]+$", dt$position)) ||
        any(!nzchar(trimws(dt$ref))) || any(!nzchar(trimws(dt$alt)))) {
      stop(modality, " SuSiE TSV has missing or invalid credible-set fields.", call. = FALSE)
    }
    dt[, modality := modality]
    data.table::setorder(dt, molecular_trait_id, cs_id, -pip, chromosome, position, ref, alt)
    dt <- dt[, .SD[1L], by = .(modality, molecular_trait_id, cs_id)]
    message(sprintf("%s: selected %d credible-set leads.", modality, nrow(dt)))
    dt
  })
  data.table::rbindlist(leads)
}

variant_match_key <- function(chromosome, position, ref, alt) {
  paste(
    sub("^chr", "", as.character(chromosome), ignore.case = TRUE),
    as.character(position),
    toupper(as.character(ref)),
    toupper(as.character(alt)),
    sep = ":"
  )
}

match_susie_leads <- function(leads, prepared) {
  required <- c("modality", "molecular_trait_id", "chromosome", "position", "ref", "alt")
  if (!all(required %in% names(leads))) {
    stop("SuSiE lead table has missing columns.", call. = FALSE)
  }
  metadata <- prepared$phenotype_metadata
  outcomes <- paste(metadata$modality, metadata$phenotype_id, sep = "::")
  gene_leads <- leads[paste(leads$modality, leads$molecular_trait_id, sep = "::") %in% outcomes, ]
  variants <- prepared$variant_metadata
  variant_keys <- variant_match_key(variants$CHROM, variants$POS, variants$REF, variants$ALT)
  if (anyDuplicated(variant_keys)) {
    stop("Prepared variants are ambiguous after chromosome normalization.", call. = FALSE)
  }
  lead_keys <- variant_match_key(
    gene_leads$chromosome, gene_leads$position, gene_leads$ref, gene_leads$alt
  )
  indices <- match(lead_keys, variant_keys)
  matched <- unique(indices[!is.na(indices)])
  message(sprintf(
    "SuSiE leads: %d credible sets, %d unmatched leads, %d unique matched variants.",
    nrow(gene_leads), sum(is.na(indices)), length(matched)
  ))
  list(
    variant_ids = colnames(prepared$X)[matched],
    n_leads = nrow(gene_leads),
    n_unmatched = sum(is.na(indices)),
    n_matched_variants = length(matched)
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 4L) {
    stop("Expected expression, splicing, protein TSV paths and output path.", call. = FALSE)
  }
  paths <- stats::setNames(args[1:3], c("expression", "splicing", "protein"))
  leads <- read_susie_leads(paths)
  data.table::fwrite(leads, args[[4L]], sep = "\t")
}
