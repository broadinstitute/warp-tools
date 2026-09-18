#!/usr/bin/env Rscript
# 
# Identifies variants that overlap other variants. If overlapping variants are found,
# only the first longest allele is kept as this is likely to have the largest impact
# on read mapping.
#
# USAGE: check_variant_bounds.R S_OUT MTVCF, NUCVCF, FORCE_NO_DUPES
#
args = commandArgs(trailingOnly=TRUE)
library(data.table)

get_ranges <- function(table) {
  # Gets positions of each alleles start and end.
  table[['pos_A1_start']] <- table[['POS']]
  table[['pos_A1_end']] <- table[['POS']] + nchar(table[['REF']]) - 1
  table[['pos_A2_start']] <- table[['POS']]
  table[['pos_A2_end']] <- table[['POS']] + nchar(table[['ALT']]) - 1
  table <- table[order(table[['pos_A1_start']]),c('CHROM','pos_A1_start', 
                                                  'pos_A1_end', 'pos_A2_start', 'pos_A2_end', 'idx')]
  return(table)
}

read_vcf <- function(path) {
  # Reads a VCF. Ignores the header.
  dt <- fread(path, skip = "#CHROM", header=T)
  setnames(dt, '#CHROM','CHROM')
  dt <- as.data.frame(dt)
  dt[['idx']] <- seq_len(nrow(dt))
  return(dt)
}

write_vcf <- function(x, path, s) {
  # Writes a custom VCF.
  x_out <- x[,c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',s)]
  if(nrow(x_out) > 0) {
    x_out[['ID']] <- '.'
    x_out[['QUAL']] <- '.'
    x_out[['FILTER']] <- 'PASS'
    x_out[['INFO']] <- '.'
    x_out[['FORMAT']] <- '.'
    x_out[[s]] <- '.'
  }
  names(x_out)[1] <- '#CHROM'
  write('##fileformat=VCFv4.2', path)
  write.table(x_out, file=path, append=T, sep='\t', quote=F, row.names=F)
}

check_overlap <- function(table, make_bins=F) {
  # For each site in a given chromosome, checks if there are overlaps.
  # ONLY WORKS PER-CHROMSOME. Makes no attempt at checking if multiple
  # chromosomes are provided.
  # make_bins=T makes this function produce "overlap IDs", thus allowing the output
  # to be used to pull out all members of an overlapping set, rather than
  # just outputting the variants that overlap with some independent set.
  #
  # OUTPUT: a vector of length nrow(table) containing nonzero overlap IDs for
  # any overlapping variant sites
  if(make_bins) {
    overlaps <- rep(0, nrow(table))
  } else {
    overlaps <- rep(F, nrow(table))
  }
  if(nrow(table) < 1) return(overlaps)
  cache_idx <- 1
  bin_idx <- 1
  for(idx in seq_along(table[['pos_A1_start']])) {
    if(make_bins) {
      if(overlaps[idx] != 0) {
        this_bin_idx <- overlaps[idx]
      } else {
        this_bin_idx <- bin_idx
      }
    }
    this_start <- table[['pos_A1_start']][idx]
    this_end <- table[['pos_A1_end']][idx]
    for(idx2 in cache_idx:nrow(table)) {
      if(idx == idx2) next
      else {
        q_start <- table[['pos_A1_start']][idx2]
        q_end <- table[['pos_A1_end']][idx2]
        if(q_start > this_end) break
        else if(q_end < this_start) {
          cache_idx <- cache_idx + 1; next
        } else {
          if(make_bins) {
            overlaps[idx] <- this_bin_idx
            overlaps[idx2] <- this_bin_idx
          } else {
            overlaps[idx] <- T
            overlaps[idx2] <- T
          }
        }
      }
    }
    bin_idx <- bin_idx + 1
  }
  return(overlaps)
}

if(length(args) != 4) {
  stop('ERROR: Must have 4 input arguments in the following order: S_OUT, MTVCF, NUCVCF, FORCE_NO_DUPES')
}

mt_in <- read_vcf(args[2])
s <- names(mt_in)[10]
mt_vcf <- get_ranges(mt_in)
if(nrow(mt_vcf) > 0) {
  split_mt <- split(mt_vcf, mt_vcf[['CHROM']])
} else {
  split_mt <- list('chrM'=mt_vcf)
}
tf_overlap <- unlist(lapply(split_mt, check_overlap))
bin_id_overlap <- lapply(split_mt, check_overlap, make_bins=T)
if (sum(unlist(bin_id_overlap) != 0) != sum(tf_overlap)) {
  stop('ERROR: somehow the results from binning and T/F outputs from check_overlap disagree.')
}
if (any(tf_overlap)) { # if FORCE_NO_DUPES is enabled, fails here if there are overlaps
  if(args[4]) {
    stop(paste0('ERROR: Mito homoplasmies VCF has ', sum(tf_overlap), ' variants which are overlapping. This will produce malformed consensus results.'))
  } else {
    warning(paste0('WARNING: Mito homoplasmies VCF has ', sum(tf_overlap), ' variants which are overlapping. This will produce malformed consensus results.'))
  }
}
write.table(data.frame(x=sum(tf_overlap)), paste0(args[1], ".mtdna_consensus_overlaps.txt"), quote=F, col.names=F, row.names=F)

tabs_for_removal <- lapply(names(split_mt), function(x) { 
  # For each chromosome, for each duplicate group, flags the non-first
  # variant for removal. Ordering is based on the input VCF order.
  tab <- split_mt[[x]]
  tab[['binid']] <- bin_id_overlap[[x]]
  tab <- tab[tab[['binid']] != 0,]
  if(nrow(tab) == 0){
    return(tab)
  } else {
    split_tab <- split(tab, tab[['binid']])
    split_tab_f <- lapply(split_tab, function(x)x[x[['idx']] != min(x[['idx']]),])
    return(do.call(rbind, split_tab_f))
  }
})
tab_removal <- do.call(rbind, tabs_for_removal)
t2 <- tab_removal[,c('pos_A1_start','idx','binid')]
variants_for_removal <- merge(x=mt_in,y=t2,by="idx",all.y=TRUE)
if(nrow(t2) != nrow(variants_for_removal)) {
  stop('ERROR: Some variants are missing or duplicate in mito VCF or in VCF for removal.')
}
if(length(unique(variants_for_removal[['idx']])) != nrow(variants_for_removal)) {
  stop('ERROR: Some variant IDs were duplicated.')
}
if(any(variants_for_removal[['POS']] != variants_for_removal[['pos_A1_start']])) {
  stop('ERROR: duplicate variants do not have the same positions as the original variants.')
}
allbins <- unique(unlist(bin_id_overlap))
allbins <- allbins[allbins != 0]
if(any(!allbins %in% variants_for_removal[['binid']])) {
  stop('ERROR: all duplicate bins must be represented in removal.')
}
if(sum(tf_overlap) - length(allbins) != nrow(variants_for_removal)) {
  stop('ERROR: the number of variants to remove must be the total number of dupes - the number of bins.')
}
write_vcf(variants_for_removal, 'overlapping_variants_to_remove.vcf', s)

nuc_vcf <- get_ranges(read_vcf(args[3]))
tf_overlap <- unlist(lapply(split(nuc_vcf, nuc_vcf[['CHROM']]), check_overlap))
if (any(tf_overlap)) {
  warning(paste0('WARNING: nucDNA homoplasmies VCF has ', sum(tf_overlap), ' variants which are overlapping. This will produce malformed consensus results.'))
}
write.table(data.frame(x=sum(tf_overlap)), paste0(args[1], ".nucdna_consensus_overlaps.txt"), quote=F, col.names=F, row.names=F)
