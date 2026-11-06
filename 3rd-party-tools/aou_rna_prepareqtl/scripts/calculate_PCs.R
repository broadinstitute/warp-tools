library(tidyverse)
library(data.table)
library(PCAtools)
library(janitor)

########## PARSE COMMAND LINE ARGUMENTS ##########

option_list <- list(
    optparse::make_option(c("--bed_file"), type="character", default=NULL,
                        help="Sample to be used in processing expression marix", metavar = "type"),
    optparse::make_option(c("--output_prefix"), type="character", default=NULL,
                        help="Sample to be used in processing expression marix", metavar = "type")
  #  optparse::make_option(c("--genetic_covariates"), type="character", default=NULL,
                        #help="Sample to be used in processing expression marix", metavar = "type")

    )
opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))



genetic_PCs <- opt$genetic_covariates
bed_file <- opt$bed_file
prefix <- opt$output_prefix

phenotype_pcs_out <- paste0(prefix,'_phenotype_PCs.tsv')
#QTL_covariates <- paste0(prefix,'_QTL_covariates.tsv')
message(paste0('Writing to: ' ,phenotype_pcs_out))

######## FUNCTIONS ########
compute_pcs <- function(expression_df){

subsetted_expression_dat <- expression_df %>% select(-c(1,2,3,4))
pca_standardized <- PCAtools::pca(subsetted_expression_dat)
n_pcs <- chooseGavishDonoho( subsetted_expression_dat ,  var.explained = pca_standardized$sdev^2, noise = 1)
message(paste0('Using' , n_pcs,' PCs'))
pca_out <- pca_standardized$rotated %>% 
   data.frame() %>%
   select(1:n_pcs) %>% 
   rownames_to_column('ID') %>% 
   mutate(ID = str_remove(ID,'X'))

pca_out
}




####### ANALYSIS BEGIN ########

bed_df <-  readr::read_tsv(bed_file)
message('Bed file loaded')

message('Computing PCs')
PCA_data <- compute_pcs(bed_df)
#genetic_PCs <- fread(genetic_PCs) %>% dplyr::rename( 'ID'= 'sample_id')

message('Writing PCs')
# writes phenotype PC to output
PCA_data %>% write_tsv(phenotype_pcs_out)

# merges genetic PCs and phenotype PCs
#merged_data <- genetic_PCs %>% 
    #inner_join(PCA_data,by = 'ID') %>% 
    #dplyr::select(ID,everything()) %>% 
    #distinct() 

## formats covariates data for tensorQTL
#output_data <- merged_data %>% 
    #arrange(ID) %>% 
    #t() %>% 
    #data.frame()
    #janitor::row_to_names(row_number = 1) %>% 
    #rownames_to_column('ID')

#output_data %>% write_tsv(QTL_covariates)

