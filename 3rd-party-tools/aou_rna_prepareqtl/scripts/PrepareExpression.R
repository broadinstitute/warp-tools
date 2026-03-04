library(tidyverse)
library(data.table)
library(magrittr)
library(biomaRt)
library(optparse)
library(data.table)
library(rtracklayer)
library(RNOmni)
library(edgeR)


# use rtracklayer to import GTF file and extract TSS locations. 
# This should run on the collapsed GTF that has been generated 
# but might work with any GTF
extract_TSS_pos <- function(gencode_file){
message('Extracting TSS locations')
message('Loading GTF')
gencode_GTF <- rtracklayer::import(gencode_file) %>% data.frame()

message('Extracting GTF')
# map TSS locations based on strand
TSS_locations <- gencode_GTF  %>% 
    filter(type == 'gene'  ) %>%
    mutate(TSS = case_when(strand == '+' ~ start,TRUE ~ end)) %>% 
    dplyr::select(gene_id,TSS,seqnames) %>% 
    mutate(start = TSS -1,end = TSS) %>% 
    dplyr::select(gene_id,start,end,seqnames)
TSS_locations
}



######## COMMAND LINE ARGUMENTS ############

option_list <- list(
    optparse::make_option(c("--CountGCT"), type="character", default=NULL,
                        help="Parquet or TSV of normalzied protein expression data", metavar = "type"),
 #   optparse::make_option(c("--TPMGCT"), type="character", default=NULL,
                        #help="Parquet or TSV of normalzied protein expression data", metavar = "type"),
    optparse::make_option(c("--OutputPrefix"), type="character", default=NULL,
                        help="Prefix for output data", metavar = "type"),
    optparse::make_option(c("--AnnotationGTF"), type="character", default=NULL,
                        help="GTF file used to TSS locations for each gene", metavar = "type"),
    optparse::make_option(c("--SampleList"), type="character", default=NULL,
                        help="File containing list of samples to run processing on", metavar = "type")
)

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))


########### LOAD DATA #####################

message('Loading count data')
CountData <-  fread(opt$CountGCT,skip  =2 ,header = TRUE)

#TPMData <- fread(opt$TPMGCT,skip  =2 ,header = TRUE)

PositionTSS <- extract_TSS_pos(opt$AnnotationGTF)

SampleList <- fread(opt$SampleList,header = FALSE) %>% dplyr::rename('ID' = 1) %>% pull(ID)
nSamples <- SampleList %>% length()
message(paste0('Number of sample in sample list:',nSamples))


OutputFile <- paste0(opt$OutputPrefix,'.expression.bed.gz')
message(paste0('Writing to file: ',OutputFile ))
############# PROCESS DATA ###########
# transpose read count data such that 
# genes are column and rows are samples 
message('Transposing data')
CountDataTransposed <- CountData %>%
    dplyr::select(-Description) %>% 
    column_to_rownames('Name') %>%
    dplyr::select(any_of(SampleList)) %>% 
    t() %>% 
    data.frame()

# filter to genes where the count is greater than 6
# in atleast 20% of samples
message('Filtering expression by counts')
CountDataFiltered <- CountDataTransposed %>%
        dplyr::select(where(~ mean(.x > 6) >= 0.2)) %>% 
        t() %>% 
        data.frame()

message('Performing edgeR TMM normalization')
# Convert filtered count data to DGE list for normalization 
DataEdgeR <- edgeR::DGEList(CountDataFiltered)
DataEdgeR <- edgeR::calcNormFactors(DataEdgeR)

message('Computing CPMs')
DataCPM <- edgeR::cpm(DataEdgeR, log=FALSE) %>% data.frame() 

message('Normalizing CPMs')
NormalizedCPMs <- DataCPM %>% 
                t() %>% 
                data.frame() %>% 
                mutate(across(everything(),~RankNorm(.))) %>% 
                t() %>% 
                data.frame() %>% 
                dplyr::rename_with(~str_remove(.,'^X')) %>% 
                rownames_to_column('gene_id')

LengthNormalziedCPMS <- NormalizedCPMs %>% nrow
message(paste0('Number of genes found: ',LengthNormalziedCPMS))

message('Merging with quantifications with TSS locations')
BedNormalizedCPMs <- PositionTSS %>% 
            inner_join(NormalizedCPMs,by = 'gene_id') %>% 
            dplyr::select(seqnames,start,end,gene_id,everything()) 


LengthBedCPMs <- BedNormalizedCPMs %>% nrow
message(paste0('Number of genes found after merge: ',LengthBedCPMs))

LostGenes <- LengthNormalziedCPMS  - LengthBedCPMs
message(paste0('Gene lost after merging: ',LostGenes))

if (LostGenes > 0){
message('Please check GENCODE version since genes are lost in merging process')

}

BedNormalizedCPMs %>% arrange(seqnames,start) %>%  dplyr::rename('#chr' = 'seqnames') %>% fwrite(OutputFile,sep='\t') 




