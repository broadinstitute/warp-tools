library(tidyverse)
library(data.table)
library(optparse)
library(janitor)

######## COMMAND LINE ARGUMENTS ############
message('MergeCovariates starting')
option_list <- list(
    optparse::make_option(c("--GenotypePCs"), type="character", default=NULL,
                        help="", metavar = "type"),
    optparse::make_option(c("--MolecularPCs"), type="character", default=NULL,
                        help="", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type")

   )

message('Parsing command line arguments')
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


OutputFile <- paste0(opt$OutputPrefix,'_QTL_covariates.tsv')

message(paste0('Writing to outputfile',OutputFile))

########## LOAD DATA #################

message('Reading genetic PCs')
geneticPCs <- fread(opt$GenotypePCs)

message('Reading molecular PCs')
molecularPCs <- fread(opt$MolecularPCs)

message('Merging data')
mergedPCs <- geneticPCs %>% 
        inner_join(molecularPCs,by = c('sample_id' = 'ID')) %>% 
        dplyr::select(sample_id,everything()) %>% 
        distinct()

message('Writing to output')
mergedPCs %>% 
    arrange(sample_id) %>% 
    t() %>% 
    data.frame() %>% 
    janitor::row_to_names(row_number = 1) %>% 
    rownames_to_column('ID')  %>% 
    write_tsv(OutputFile)

