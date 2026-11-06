#!/usr/bin/env Rscript
# This script maps the header of FASTA segments to the intervals provided.
#
# Applying ExtractSequences produces an output FASTA file which contains headers
# corresponding to the target name in the applied interval file. bcftools consensus, 
# however, expects the FASTA header to have format >chr:from-to so it knows where in the
# genome the FASTA file maps to relative to the provided VCF. This script
# allows for mapping back and forth between a FASTA with the >chr:from-to format and a
# FASTA with headers corresponding to the interval names.
#
# USAGE: compatibilify_fa_intervals_consensus.R FASTA INTERVAL_LIST REVERSE(bool) OUTPUT MITOMODE
#
# Typically run once after ExtractSequences to prepare for bcftools consensus, and then again
# in reverse mode to map the consensus segments back to interval names.
# Throws an error if the interval is not found.
#
# 221104: updated to work with empty interval lists and empty fasta files
#
args = commandArgs(trailingOnly=TRUE)

generate_map <- function(interval, map_back) {
  # Generates a map between interval names and chr:from-to formatted names.
  # If map_back is enabled, the outputted vector maps from chr:from-to to names.
  spl_interval <- strsplit(interval,'\t')
  nms <- sapply(spl_interval, function(x)x[5])
  target <- sapply(spl_interval, function(x)paste0(x[1],':',x[2],'-',x[3]))
  if(any(duplicated(nms))) stop('ERROR: Duplicate interval names not allowed.')
  if(any(duplicated(target))) stop('ERROR: Duplicated ranges for the interval not allowed.')
  if(map_back) {
    names(nms) <- target
    return(nms)
  } else {
    names(target) <- nms
    return(target)
  }

}

rename_fasta <- function(fasta, name_map, require_all) {
  # After validating that all fasta segments map each to exactly one interval
  # file segment and vice versa, renames the fasta according to name_map.
  fasta_out <- fasta
  to_change <- grep('^>', fasta_out)
  to_change_conv <- gsub('^>','',fasta_out[to_change])
  if(!all(to_change_conv %in% names(name_map))) {
    stop('ERROR: All FASTA segments must be found in interval file.')
  }
  if(require_all) {
    if(!all(names(name_map) %in% to_change_conv)) {
      stop('ERROR: All interval file segments must be found in FASTA file.')
    }
  }
  fasta_out[to_change] <- paste0('>', name_map[to_change_conv])
  return(fasta_out)
}

if(length(args) != 5) {
  stop('ERROR: Must have 5 input arguments in the following order: FASTA, INTERVAL_LIST, REVERSE(bool), OUTPUT')
}

fast <- readLines(args[1])
interval <- readLines(args[2])
interval <- interval[grep('^@',interval, invert=T)]

if(length(fast) == 0) {
  print('FASTA is empty.')
  if(length(interval) != 0) {
    stop('ERROR: cannot have an empty fasta with a nonempty interval.')
  }
  if(args[5]) {
    stop('ERROR: cannot have an empty fasta with mito mode enabled.')
  }
  writeLines(fast, args[4])
} else {
  if(args[5]) {
    # If MITOMODE is enabled, validates the provided interval.
    tf_chr <- sapply(strsplit(interval, '\t'), function(x)x[1] == 'chrM')
    tf_name <- sapply(strsplit(interval, '\t'), function(x)x[5] == 'chrM')
    tf_start <- sapply(strsplit(interval, '\t'), function(x)x[2] == '1')
    if(!((length(interval) == 1) & all(tf_chr) & all(tf_start) & all(tf_name))) {
      stop('ERROR: If mito mode is enabled for interals, there must be a single interval starting at 1 with chrM.')
    }
  }
  name_map <- generate_map(interval, args[3])
  renamed_fasta <- rename_fasta(fast, name_map, T)
  writeLines(renamed_fasta, args[4])
}