#!/usr/bin/env Rscript

# Prepare covariates for AoU proteomics analysis from Olink NPX data.
# This script converts the notebook workflow into a parameterized CLI suitable for WDL 1.0.

suppressPackageStartupMessages({
	library(data.table)
	library(dplyr)
	library(tidyr)
	library(tibble)
	library(stringr)
	library(WGCNA)
	library(biomaRt)
	library(rtracklayer)
	library(PCAtools)
	library(RNOmni)
	library(janitor)
})

# ------------------------------
# CLI parsing and validation
# ------------------------------

parse_args <- function(args) {
	out <- list()
	i <- 1
	while (i <= length(args)) {
		arg <- args[[i]]
		if (!startsWith(arg, "--")) {
			stop("Unexpected argument format: ", arg)
		}

		keyval <- sub("^--", "", arg)
		if (grepl("=", keyval, fixed = TRUE)) {
			parts <- strsplit(keyval, "=", fixed = TRUE)[[1]]
			key <- parts[[1]]
			val <- paste(parts[-1], collapse = "=")
			out[[key]] <- val
			i <- i + 1
		} else {
			key <- keyval
			if (i == length(args) || startsWith(args[[i + 1]], "--")) {
				out[[key]] <- TRUE
				i <- i + 1
			} else {
				out[[key]] <- args[[i + 1]]
				i <- i + 2
			}
		}
	}
	out
}

required_args <- c(
	"comb_psam",
	"olink_normalized",
	"ancestry_preds",
	"gencode_gtf",
	"genetic_pc_source_path",
	"genetic_pc_basename",
	"bed_bucket_path",
	"phenotype_pc_bucket_path",
	"covariate_bucket_path"
)

args <- parse_args(commandArgs(trailingOnly = TRUE))

missing_args <- required_args[!vapply(required_args, function(x) !is.null(args[[x]]), logical(1))]
if (length(missing_args) > 0) {
	stop("Missing required args: ", paste(missing_args, collapse = ", "))
}

output_dir <- if (!is.null(args$output_dir)) args$output_dir else "."
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
setwd(output_dir)

output_prefix <- if (!is.null(args$output_prefix)) args$output_prefix else "AoU"
phenotype_suffix <- if (!is.null(args$phenotype_suffix)) args$phenotype_suffix else "10k_pQTL_phenotype_PCs.tsv"
bed_suffix <- if (!is.null(args$bed_suffix)) args$bed_suffix else "10k_pQTL.bed.gz"
covar_sample_label <- if (!is.null(args$covar_sample_label)) args$covar_sample_label else "10750"
connectivity_basename <- if (!is.null(args$connectivity_basename)) args$connectivity_basename else "olink_10k"

# ------------------------------
# Utility helpers
# ------------------------------

# Localize a GCS object if needed, otherwise return the existing local path.
localize_path <- function(path) {
	if (startsWith(path, "gs://")) {
		local_name <- basename(path)
		cmd <- sprintf("gsutil -m cp %s %s", shQuote(path), shQuote(local_name))
		status <- system(cmd)
		if (status != 0) stop("Failed to localize: ", path)
		return(local_name)
	}
	if (!file.exists(path)) stop("Local input does not exist: ", path)
	path
}

# Copy a local output file to a target bucket path.
copy_to_bucket <- function(local_file, bucket_path) {
	if (is.null(bucket_path) || bucket_path == "") return(invisible(NULL))
	cmd <- sprintf("gsutil cp %s %s", shQuote(local_file), shQuote(bucket_path))
	status <- system(cmd)
	if (status != 0) stop("Failed to copy to bucket: ", local_file, " -> ", bucket_path)
}

# Build a file path from a directory/bucket prefix and filename.
build_path <- function(prefix, filename) {
	if (is.null(prefix) || prefix == "") return(filename)
	file.path(prefix, filename)
}

# ------------------------------
# Input localization
# ------------------------------

comb_psam <- localize_path(args$comb_psam)
olink_normalized <- localize_path(args$olink_normalized)
ancestry_preds <- localize_path(args$ancestry_preds)
gencode_gtf <- localize_path(args$gencode_gtf)

# Build ancestry genetic-PC paths from a single source path and a basename suffix.
# Example: genetic_pc_source_path=gs://bucket/path and genetic_pc_basename=genetic_PCs.tsv
# yields gs://bucket/path/MID_genetic_PCs.tsv, AFR_genetic_PCs.tsv, etc.
genetic_pc_source_path <- args$genetic_pc_source_path
genetic_pc_basename <- args$genetic_pc_basename

genetic_pc_input_paths <- list(
	comb = build_path(genetic_pc_source_path, paste0("COMB_", genetic_pc_basename)),
	mid = build_path(genetic_pc_source_path, paste0("MID_", genetic_pc_basename)),
	afr = build_path(genetic_pc_source_path, paste0("AFR_", genetic_pc_basename)),
	eas = build_path(genetic_pc_source_path, paste0("EAS_", genetic_pc_basename)),
	eur = build_path(genetic_pc_source_path, paste0("EUR_", genetic_pc_basename)),
	sas = build_path(genetic_pc_source_path, paste0("SAS_", genetic_pc_basename)),
	amr = build_path(genetic_pc_source_path, paste0("AMR_", genetic_pc_basename))
)

genetic_pc_local_paths <- lapply(genetic_pc_input_paths, localize_path)

# ------------------------------
# Load primary proteomics inputs
# ------------------------------

multiomics_sample_list <- fread(comb_psam) %>% dplyr::rename("ResearchID" = 1)
olink_df <- fread(olink_normalized) %>%
	filter(ResearchID %in% multiomics_sample_list$ResearchID)

# ------------------------------
# Map UniProt IDs to Ensembl gene IDs
# ------------------------------

protein_list <- olink_df %>%
	dplyr::select(UniProt) %>%
	distinct() %>%
	pull(UniProt)

mart <- useEnsembl(
	biomart = "ensembl",
	dataset = "hsapiens_gene_ensembl",
	host = "https://www.ensembl.org"
)

conversion_list <- getBM(
	attributes = c("ensembl_gene_id_version", "uniprot_gn_id"),
	filters = "uniprot_gn_id",
	values = protein_list,
	mart = mart
) %>%
	dplyr::rename(
		gene_id = ensembl_gene_id_version,
		UniProt = uniprot_gn_id
	)

# ------------------------------
# Compute assay and per-person missingness
# ------------------------------

missingness_per_assay <- olink_df %>%
	group_by(Assay) %>%
	summarize(n = 100 * sum(is.na(PCNormalizedNPX)) / dplyr::n(), .groups = "drop")

missingness_per_person <- olink_df %>%
	group_by(SampleID) %>%
	summarize(n = sum(is.na(PCNormalizedNPX)), .groups = "drop")

person_keep_list <- missingness_per_person %>% filter(n < 1)
message("Number of individuals with no missing assays: ", nrow(person_keep_list))

# ------------------------------
# Convert to wide matrix for network outlier detection
# ------------------------------

NPX_data_wide <- olink_df %>%
	filter(SampleID %in% person_keep_list$SampleID) %>%
	group_by(ResearchID, OlinkID) %>%
	filter(row_number() == 1) %>%
	ungroup() %>%
	filter(Assay != "GBP1" & Assay != "MAP2K1") %>%
	dplyr::select(ResearchID, NPX, UniProt) %>%
	pivot_wider(names_from = UniProt, values_from = NPX) %>%
	column_to_rownames("ResearchID")

# ------------------------------
# Compute connectivity Z-scores and identify outlier samples
# ------------------------------

norm_adj <- 0.5 + 0.5 * bicor(NPX_data_wide %>% t())
net_summary <- fundamentalNetworkConcepts(norm_adj)
net_connectivity <- net_summary$Connectivity

connectivity_zscore <- ((net_connectivity - mean(net_connectivity)) / sd(net_connectivity)) %>%
	data.frame() %>%
	dplyr::rename("Z_score" = 1) %>%
	rownames_to_column("ResearchID")

connectivity_zscore_outliers <- connectivity_zscore %>% filter(Z_score < -3)
message("Number of connectivity Z-score outliers: ", nrow(connectivity_zscore_outliers))

zscore_file <- paste0(connectivity_basename, "_connectivity_z_scores.tsv")
outlier_file <- paste0(connectivity_basename, "_connectivity_z_score_outliers.tsv")
sample_list_file <- paste0(connectivity_basename, "_sample_list.tsv")

connectivity_zscore %>% write_tsv(zscore_file)
connectivity_zscore_outliers %>% dplyr::select(ResearchID) %>%
	write_tsv(outlier_file, col_names = FALSE)

copy_to_bucket(outlier_file, args$covariate_bucket_path)
copy_to_bucket(zscore_file, args$covariate_bucket_path)

# ------------------------------
# Remove outlier samples and write the retained sample list
# ------------------------------

connectivity_zscore_outliers <- fread(
	outlier_file,
	header = FALSE,
	col.names = "ResearchID"
)

PCNPX_wide_filtered <- NPX_data_wide %>%
	rownames_to_column("ResearchID") %>%
	filter(!ResearchID %in% connectivity_zscore_outliers$ResearchID) %>%
	column_to_rownames("ResearchID")

message("Number of samples after removing outliers: ", nrow(PCNPX_wide_filtered))
rownames(PCNPX_wide_filtered) %>%
	data.frame() %>%
	write_tsv(sample_list_file, col_names = FALSE)

# ------------------------------
# Evaluate optional PASS-only filtering impact (kept for parity with notebook)
# ------------------------------

df_qc <- olink_df %>%
	filter(SampleID %in% person_keep_list$SampleID)

n_samples_before <- df_qc %>% distinct(SampleID) %>% nrow()
n_assays_before <- df_qc %>% distinct(OlinkID) %>% nrow()
n_rows_before <- nrow(df_qc)

df_pass <- df_qc %>%
	filter(AssayQC == "PASS", SampleQC == "PASS")

n_samples_after <- df_pass %>% distinct(SampleID) %>% nrow()
n_assays_after <- df_pass %>% distinct(OlinkID) %>% nrow()
n_rows_after <- nrow(df_pass)

message(
	"PASS/PASS filter impact - samples lost: ", n_samples_before - n_samples_after,
	", assays lost: ", n_assays_before - n_assays_after,
	", measurements lost: ", n_rows_before - n_rows_after
)

# ------------------------------
# Build gene-level matrix with TSS coordinates from GTF
# ------------------------------

# Remap UniProt IDs to Ensembl IDs using Ensembl version 113 (matches notebook).
ensembl <- useEnsembl(
	biomart = "genes",
	dataset = "hsapiens_gene_ensembl",
	version = 113
)

conversion_list <- biomaRt::getBM(
	attributes = c("ensembl_gene_id_version", "uniprot_gn_id"),
	filters = "uniprot_gn_id",
	values = protein_list,
	mart = ensembl
) %>%
	dplyr::rename(
		gene_id = ensembl_gene_id_version,
		UniProt = uniprot_gn_id
	)

gencode_GTF <- rtracklayer::import(gencode_gtf) %>% data.frame()

TSS_locations <- gencode_GTF %>%
	filter(type == "gene") %>%
	mutate(TSS = case_when(strand == "+" ~ start, TRUE ~ end)) %>%
	dplyr::select(gene_id, TSS, seqnames) %>%
	mutate(start = TSS - 1, end = TSS) %>%
	dplyr::select(gene_id, start, end, seqnames)

gene_id_merged_NPX_data <- PCNPX_wide_filtered %>%
	t() %>%
	data.frame() %>%
	rownames_to_column("UniProt") %>%
	left_join(conversion_list, by = "UniProt") %>%
	left_join(TSS_locations, by = "gene_id") %>%
	dplyr::select(seqnames, start, end, UniProt, gene_id, everything()) %>%
	filter(!is.na(seqnames)) %>%
	dplyr::rename_with(~str_remove(., "X")) %>%
	mutate(gene_id = paste0(UniProt, "_", gene_id)) %>%
	dplyr::select(-UniProt) %>%
	arrange(seqnames, start)

# ------------------------------
# Load ancestry predictions
# ------------------------------

ancestry_df <- fread(ancestry_preds) %>%
	dplyr::select(research_id, ancestry_pred_other)

# ------------------------------
# Generate per-population phenotype bed and phenotype PCs
# ------------------------------

write_covars_to_bucket <- function(
	ancestry_df,
	group_label,
	bed_data,
	bed_bucket_path = "",
	covar_bucket_path = "",
	covariate_suffix = "",
	bed_suffix = "",
	filter_group = TRUE,
	output_prefix = "AoU"
) {
	message("Running on population: ", group_label)

	bed_file_name <- paste0(output_prefix, "_", group_label, "_", bed_suffix)
	covar_file_name <- paste0(output_prefix, "_", group_label, "_", covariate_suffix)

	if (filter_group) {
		group_ids <- ancestry_df %>%
			filter(ancestry_pred_other == group_label) %>%
			mutate(research_id = as.character(research_id)) %>%
			pull(research_id)
	} else {
		group_ids <- ancestry_df %>%
			mutate(research_id = as.character(research_id)) %>%
			pull(research_id)
	}

	subset_bed <- bed_data %>% dplyr::select(1:4, any_of(group_ids))
	meta_data <- subset_bed %>% dplyr::select(1:4)

	normalized_quantifications <- subset_bed %>%
		dplyr::select(any_of(group_ids)) %>%
		t() %>%
		data.frame() %>%
		mutate(across(everything(), ~RNOmni::RankNorm(.))) %>%
		t() %>%
		data.frame() %>%
		arrange()

	rownames(normalized_quantifications) <- NULL

	merged_data <- bind_cols(meta_data, normalized_quantifications) %>%
		dplyr::rename_with(~str_remove(., "X")) %>%
		dplyr::rename("#chr" = "seqnames")

	merged_data %>% fwrite(bed_file_name, sep = "\t")

	PCA_res <- merged_data %>%
		data.frame() %>%
		column_to_rownames("gene_id") %>%
		dplyr::select(-1, -2, -3) %>%
		PCAtools::pca()

	num_PCs <- PCAtools::chooseGavishDonoho(
		normalized_quantifications,
		var.explained = PCA_res$sdev^2,
		noise = 1
	)

	message("Number of PCs chosen for ", group_label, ": ", num_PCs)

	PCA_calculations <- PCA_res$rotated %>%
		dplyr::select(1:num_PCs) %>%
		data.frame() %>%
		rownames_to_column("research_id") %>%
		mutate(research_id = str_remove(research_id, "X"))

	PCA_calculations %>% fwrite(covar_file_name, sep = "\t")

	copy_to_bucket(covar_file_name, covar_bucket_path)
	copy_to_bucket(bed_file_name, bed_bucket_path)

	number_samples_bed <- subset_bed %>% dplyr::select(-1, -2, -3, -4) %>% ncol()
	message("Number of samples in subset bed for ", group_label, ": ", number_samples_bed)

	list(bed_file = bed_file_name, covar_file = covar_file_name)
}

populations <- c("mid", "afr", "eur", "sas", "eas")
generated_files <- list()

for (pop in populations) {
	generated_files[[pop]] <- write_covars_to_bucket(
		ancestry_df = ancestry_df,
		group_label = pop,
		bed_data = gene_id_merged_NPX_data,
		bed_bucket_path = args$bed_bucket_path,
		covar_bucket_path = args$phenotype_pc_bucket_path,
		covariate_suffix = phenotype_suffix,
		bed_suffix = bed_suffix,
		filter_group = TRUE,
		output_prefix = output_prefix
	)
}

generated_files[["amr"]] <- write_covars_to_bucket(
	ancestry_df = ancestry_df,
	group_label = "amr",
	bed_data = gene_id_merged_NPX_data %>% distinct(),
	bed_bucket_path = args$bed_bucket_path,
	covar_bucket_path = args$phenotype_pc_bucket_path,
	covariate_suffix = phenotype_suffix,
	bed_suffix = bed_suffix,
	filter_group = TRUE,
	output_prefix = output_prefix
)

generated_files[["comb"]] <- write_covars_to_bucket(
	ancestry_df = ancestry_df,
	group_label = "comb",
	bed_data = gene_id_merged_NPX_data,
	bed_bucket_path = args$bed_bucket_path,
	covar_bucket_path = args$phenotype_pc_bucket_path,
	covariate_suffix = phenotype_suffix,
	bed_suffix = bed_suffix,
	filter_group = FALSE,
	output_prefix = output_prefix
)

# ------------------------------
# Merge phenotype PCs with genetic PCs per population
# ------------------------------

comb_genetic_PCs <- fread(genetic_pc_local_paths$comb)
mid_genetic_PCs <- fread(genetic_pc_local_paths$mid)
afr_genetic_PCs <- fread(genetic_pc_local_paths$afr)
eas_genetic_PCs <- fread(genetic_pc_local_paths$eas)
sas_genetic_PCs <- fread(genetic_pc_local_paths$sas)
eur_genetic_PCs <- fread(genetic_pc_local_paths$eur)
amr_genetic_PCs <- fread(genetic_pc_local_paths$amr)

comb_phenotype_PCs <- fread(generated_files[["comb"]]$covar_file)
mid_phenotype_PCs <- fread(generated_files[["mid"]]$covar_file)
eur_phenotype_PCs <- fread(generated_files[["eur"]]$covar_file)
sas_phenotype_PCs <- fread(generated_files[["sas"]]$covar_file)
afr_phenotype_PCs <- fread(generated_files[["afr"]]$covar_file)
eas_phenotype_PCs <- fread(generated_files[["eas"]]$covar_file)
amr_phenotype_PCs <- fread(generated_files[["amr"]]$covar_file)

comb_merged_covariates <- comb_genetic_PCs %>%
	inner_join(comb_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

mid_merged_covariates <- mid_genetic_PCs %>%
	inner_join(mid_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

sas_merged_covariates <- sas_genetic_PCs %>%
	inner_join(sas_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

afr_merged_covariates <- afr_genetic_PCs %>%
	inner_join(afr_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

eur_merged_covariates <- eur_genetic_PCs %>%
	inner_join(eur_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

eas_merged_covariates <- eas_genetic_PCs %>%
	inner_join(eas_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

amr_merged_covariates <- amr_genetic_PCs %>%
	inner_join(amr_phenotype_PCs, by = c("sample_id" = "research_id")) %>%
	dplyr::select(sample_id, everything()) %>%
	distinct()

# ------------------------------
# Write final transposed covariate outputs and copy to bucket
# ------------------------------

write_final_covar <- function(df, pop) {
	out_name <- paste0(output_prefix, "_", toupper(pop), "_", covar_sample_label, "_pQTL_covars.tsv")
	df %>%
		arrange(sample_id) %>%
		t() %>%
		data.frame() %>%
		janitor::row_to_names(row_number = 1) %>%
		rownames_to_column("ID") %>%
		write_tsv(out_name)
	copy_to_bucket(out_name, args$covariate_bucket_path)
	out_name
}

final_files <- c(
	write_final_covar(amr_merged_covariates, "amr"),
	write_final_covar(comb_merged_covariates, "comb"),
	write_final_covar(afr_merged_covariates, "afr"),
	write_final_covar(mid_merged_covariates, "mid"),
	write_final_covar(sas_merged_covariates, "sas"),
	write_final_covar(eas_merged_covariates, "eas"),
	write_final_covar(eur_merged_covariates, "eur")
)

# ------------------------------
# Emit run manifest for WDL collection/debugging
# ------------------------------

manifest <- list(
	genetic_pc_inputs = unlist(genetic_pc_input_paths),
	connectivity_zscores = zscore_file,
	connectivity_outliers = outlier_file,
	retained_sample_list = sample_list_file,
	phenotype_pc_files = vapply(generated_files, function(x) x$covar_file, character(1)),
	bed_files = vapply(generated_files, function(x) x$bed_file, character(1)),
	final_covariates = final_files
)

saveRDS(manifest, file = "prep_covariates_outputs.rds")
message("Finished prep_covariates workflow.")
