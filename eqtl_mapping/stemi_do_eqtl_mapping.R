# load libraries
library(MatrixEQTL)
library(data.table)

# functions
do_QTL_mapping <- function(
  SNP_file_name, # Genotype file name
  expression_file_name, # Gene expression file name
  snpspos, # dataframe containing the snp positions
  genepos, # dataframe containing the gene positions
  output_file_name_cis, # Output file name
  covariates_file_name=character(), # Covariates file name
  output_file_name_tra=tempfile(), # Output file name
  pvOutputThreshold_cis=1, # Only associations significant at this level will be saved
  pvOutputThreshold_tra=0, # Only associations significant at this level will be saved
  useModel=modelLINEAR, # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  cisDist=1e6, # Distance for local gene-SNP pairs
  errorCovariance=numeric(), # Error covariance matrix
  snps_fileDelimiter='\t',# the TAB character
  snps_fileOmitCharacters='NA', # denote missing values;
  snps_fileSkipRows=1, # one row of column labels
  snps_fileSkipColumns=1, # one column of row labels
  snps_fileSliceSize=2000, # read file in slices of 2,000 rows
  gene_fileDelimiter='\t',# the TAB character
  gene_fileOmitCharacters='NA', # denote missing values;
  gene_fileSkipRows=1, # one row of column labels
  gene_fileSkipColumns=1, # one column of row labels
  gene_fileSliceSize=2000, # read file in slices of 2,000 rows
  cvrt_fileDelimiter='\t',# the TAB character
  cvrt_fileOmitCharacters='NA', # denote missing values;
  cvrt_fileSkipRows=1, # one row of column labels
  cvrt_fileSkipColumns=1, # one column of row labels
  verbose=TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
){
  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = snps_fileDelimiter      # the TAB character
  snps$fileOmitCharacters = snps_fileOmitCharacters; # denote missing values;
  snps$fileSkipRows = snps_fileSkipRows;          # one row of column labels
  snps$fileSkipColumns = snps_fileSkipColumns;       # one column of row labels
  snps$fileSliceSize = snps_fileSliceSize;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = gene_fileDelimiter;      # the TAB character
  gene$fileOmitCharacters = gene_fileOmitCharacters; # denote missing values;
  gene$fileSkipRows = gene_fileSkipRows;          # one row of column labels
  gene$fileSkipColumns = gene_fileSkipColumns;       # one column of row labels
  gene$fileSliceSize = gene_fileSliceSize;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = cvrt_fileDelimiter;      # the TAB character
  cvrt$fileOmitCharacters = cvrt_fileOmitCharacters; # denote missing values;
  cvrt$fileSkipRows = cvrt_fileSkipRows;          # one row of column labels
  cvrt$fileSkipColumns = cvrt_fileSkipColumns;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  # run the actual application
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = verbose,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = min.pv.by.genesnp,
    noFDRsaveMemory = noFDRsaveMemory);
  # release resources for the temporary file
}

# not yet configurable constants
maf <- 0.1
cisDist <- 1e5

# read the command line arguments
args <- commandArgs(trailingOnly=TRUE)

# grab the snps and features
features_loc <- args[1]
snps_loc <- args[2]
# grab positions
snps_location_file_name <- args[3]
gene_location_file_name <- args[4]
# get cis mapping output loc
output_file_name_cis <- args[5]
# get covariate data
covariates_file_name <- args[6]
# TODO expand further

test=T
if(test){
  features_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_all_lowerres_20210301_metaqtl/Baseline/monocyte_expression.tsv'
  snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  output_file_name_cis<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301/Baseline/monocyte.cis.tsv'
  covariates_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_all_lowerres_20210301/Baseline/monocyte_metadata.chem.tsv'
}



# read covariate data
covariates <- fread(covariates_file_name, sep = '\t', header = T, stringsAsFactors=FALSE)
# read the expression data
expressions <- fread(features_loc, sep = '\t', header=T, stringsAsFactors=FALSE)
# read the SNPs
snps <- fread(snps_loc, sep = '\t', header = T, stringsAsFactors=FALSE)
# get the participants that we have both expression and snps data for
participants <- intersect(colnames(snps)[2:ncol(snps)], colnames(expressions)[2:ncol(expressions)])
# get also overlap with the covariates data
participants <- intersect(participants, colnames(covariates)[2:ncol(covariates)])
# perform subsetting
snps <- snps[, c('id', participants), with = F]
expressions <- expressions[, c('id', participants), with = F]
covariates <- covariates[, c('id', participants), with = F]
# filter snps by maf, which of course can work both ways
maf_reverse <- 1-maf
snps <- snps[rowSums(snps[, 2:ncol(snps)])/(ncol(snps)-1)/2 >= maf & rowSums(snps[, 2:ncol(snps)])/(ncol(snps)-1)/2 <= maf_reverse, ,]
# now write the filtered files to temporary storage
SNP_file_name <- tempfile()
expression_file_name <- tempfile()
covariates_file_name <- tempfile()
write.table(snps, SNP_file_name, sep = '\t', row.names = F, col.names = T)
write.table(expressions, expression_file_name, sep = '\t', row.names = F, col.names = T)
write.table(covariates, covariates_file_name, sep = '\t', row.names = F, col.names = T)
# next read the position files
snpspos = fread(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = fread(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
# finally call the function
do_QTL_mapping(
  SNP_file_name=SNP_file_name, # Genotype file name
  expression_file_name=expression_file_name, # Gene expression file name
  snpspos=snpspos, # dataframe containing the snp positions
  genepos=genepos, # dataframe containing the gene positions
  output_file_name_cis=output_file_name_cis, # Output file name
  covariates_file_name=covariates_file_name # Covariates file name
)



