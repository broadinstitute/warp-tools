#!/usr/bin/env Rscript
# 
# Identifies variants that should fail because they span the edge of an interval
# or because they overlap other variants. If overlapping variants are found,
# only the first longest allele is kept as this is likely to have the largest impact
# on read mapping.
#
# USAGE: check_variant_bounds.R VCF INTERVALS REJECTS
#
# 221104: updated to deal with 0-length variant callsets
#
args = commandArgs(trailingOnly=TRUE)
library(data.table)

## Helper functions
get_ranges <- function(table) {
  # Gets positions of each alleles start and end.
  table[['pos_A1_start']] <- table[['POS']]
  table[['pos_A1_end']] <- table[['POS']] + nchar(table[['REF']]) - 1
  table[['pos_A2_start']] <- table[['POS']]
  table[['pos_A2_end']] <- table[['POS']] + nchar(table[['ALT']]) - 1
  table <- table[order(table[['pos_A1_start']]),c('CHROM','REF','ALT','pos_A1_start', 
                                                  'pos_A1_end', 'pos_A2_start', 'pos_A2_end')]
  return(table)
}

read_vcf <- function(path) {
  # Reads a VCF. Saves the header and the VCF file.
  header <- readLines(path)
  header <- header[1:(grep('^#CHROM', header)[1]-1)]
  dt <- fread(path, skip = "#CHROM", header=T)
  setnames(dt, '#CHROM','CHROM')
  return(list('table' = as.data.frame(dt), 'header' = header))
}

write_vcf <- function(tab, header, path) {
  # Outputs a VCF file and includes a provided header.
  nm <- names(tab)
  nm[nm == 'CHROM'] <- '#CHROM'
  header <- c(header, paste0(nm, collapse = '\t'))
  if(nrow(tab) == 0) {
    output <- header
  } else {
    data <- sapply(seq_len(nrow(tab)), function(x)paste0(unlist(tab[x,]), collapse='\t'))
    output <- c(header, data)
  }
  writeLines(output, path)
}

assign_target <- function(chr, start, targets) {
  # Maps a site to a target table. Throws an error if the site is not mapped
  # to exactly one target.
  tab_out <- targets[((targets[['chr']] == chr) & 
                       (targets[['start']] <= start) &
                       (targets[['end']] >= start)),]
  if(nrow(tab_out) != 1) {
    stop('ERROR: should have exactly 1 match.')
  }
  return(data.frame('name' = tab_out[['name']], 
                    'start' = tab_out[['start']], 
                    'end' = tab_out[['end']]))
}

make_var_id <- function(tab, pos) {
  # Convert an imported VCF file into variant ID format.
  return(paste(tab[['CHROM']], tab[[pos]], 
               tab[['REF']], tab[['ALT']], sep=':'))
}

check_overlap <- function(table) {
  # For each site in a given chromosome, checks if there are overlaps.
  # ONLY WORKS PER-CHROMSOME. Makes no attempt at checking if multiple
  # chromosomes are provided.
  # This function always produces "overlap IDs", thus allowing the output
  # to be used to pull out all members of an overlapping set, rather than
  # just outputting the variants that overlap with some independent set.
  #
  # OUTPUT: a vector of length nrow(table) containing nonzero overlap IDs for
  # any overlapping variant sites
  overlaps <- rep(0, nrow(table))
  id <- 1
  incr <- F
  if(nrow(table) < 1) return(overlaps)
  cache_idx <- 1
  for(idx in seq_along(table[['pos_A1_start']])) {
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
          overlaps[idx] <- id
          overlaps[idx2] <- id
          incr <- T
        }
      }
    }
    if (incr) { id <- id + 1; incr <- F }
  }
  return(overlaps)
}

make_failures <- function(tab) {
  # For each set of overlapping variants, keep only one of the variants
  # with the longest reference allele. 
  # OUTPUT: failed variants as a vector in CHR:POS:REF:ALT format.
  #
  # TODO: currently does not do this per-chromosome, so theoretically
  # if two sets have failed across two chromosomes these may be considered one
  # set here. This is a conservative issue, producing loss of n-1 extra variants
  # for each n sets that are merged. In practice, this issue has not arisen. This
  # issue will not impact single-contig lists.
  print('nucDNA overlapping variants:')
  print(tab)
  if(nrow(tab) == 0) return(vector('character',0))
  else {
    lens <- sapply(split(tab, tab[['dupe_idx']]), nrow)
    if(any(lens == 1)) {
      stop('ERROR: variants entering make_failures should be in sets of at least 2.')
    }
    new_tabs <- lapply(split(tab, tab[['dupe_idx']]), function(x) {
      x[['len']] <- x[['pos_A1_end']] - x[['pos_A1_start']]
      idx_max <- which(x[['len']] == max(x[['len']]))[1]
      tf_not_max <- !(1:nrow(x) %in% idx_max)
      return(x[tf_not_max,])
    })
    lens_new <- sapply(new_tabs, nrow)
    if (!all(lens_new + 1 == lens)) {
      stop('ERROR: issues with make_failures result in unexpected number of drops.')
    }
    return(make_var_id(rbindlist(new_tabs), 'pos_A1_start'))
  }
}

## Verify arguments
if(length(args) != 3) {
  stop('ERROR: Must have 3 input arguments in the following order: VCF INTERVALS REJ')
}

## Get failures which are out of bounds of the intervals file
nuc_vcf <- read_vcf(args[1])

## If there are no variants, then skip this whole thing and output the original VCF
if (nrow(nuc_vcf[[1]]) == 0) {
  print('No variant calls found.')
  write_vcf(nuc_vcf[[1]], nuc_vcf[[2]], args[3])

} else {
  nuc_ranges <- get_ranges(nuc_vcf[[1]]) # ranges of each variant in the VCF
  nuc_intervals <- read.table(args[2], comment.char='@', sep='\t', row.names = NULL,
                              stringsAsFactors=F, col.names=c('chr','start','end','strand','name'))
  if (nrow(nuc_intervals) == 0) {
    stop('ERROR: Cannot square getting an intervals file with 0 records and actual variant calls.')
  }
  s_intervals <- split(nuc_intervals,nuc_intervals[['chr']])
  res_tab <- lapply(1:nrow(nuc_ranges), function(x) { # maps of each variant in the VCF to the interval from which it comes from
    assign_target(nuc_ranges[['CHROM']][x], nuc_ranges[['pos_A1_start']][x], 
                  s_intervals[[nuc_ranges[['CHROM']][x]]]) 
  })
  res_tab <- rbindlist(res_tab)
  failures <- nuc_ranges[(nuc_ranges[['pos_A1_start']] < res_tab[['start']]) | 
                          (nuc_ranges[['pos_A1_end']] > res_tab[['end']]),]
  failures <- make_var_id(failures, 'pos_A1_start')

  ## Get sites which are overlapping
  by_chr <- split(nuc_ranges, nuc_ranges[['CHROM']])
  idx_dupe <- lapply(by_chr, check_overlap)
  pos_dupe <- sapply(idx_dupe, function(x)any(x != 0))
  filtered_dupe_list <- rbindlist(lapply(names(pos_dupe[pos_dupe]), function(x) {
    tab <- by_chr[[x]]
    tab[['dupe_idx']] <- idx_dupe[[x]]
    return(tab)
  }))
  dupe_fails <- make_failures(filtered_dupe_list[filtered_dupe_list[['dupe_idx']] != 0,])
  failures <- unique(c(dupe_fails, failures))

  ## Output failed variants
  print('Removing the following variants:')
  print(failures) # contains duplicates and variants that span the boundary of the interval
  write_vcf(nuc_vcf[[1]][make_var_id(nuc_vcf[[1]], 'POS') %in% failures,], nuc_vcf[[2]], args[3])
}

