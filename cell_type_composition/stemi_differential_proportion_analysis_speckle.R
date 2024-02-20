#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_differential_proportion_analysis_speckle.R
# Function: perform differential proportion analysis with speckle
############################################################################################################################

####################
# libraries        #
####################

library(speckle) #BiocManager::install(c("CellBench", "BiocStyle", "scater", "org.Mm.eg.db"))
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)


####################
# Functions        #
####################

#' perform differential proportion analysis
#' 
#' @param seurat_object The Seurat object to perform differential proportion on, supply this or straight-up metadata
#' @param metadata The metadata to perform differential proportion on, supply this or the Seurat object
#' @param cell_type_column the column in the metadata describing the celltype
#' @param variable_of_interest the variable we are interested in for DPA
#' @param participant_column the column denoting the participant the cells came from
#' @param covariates vector of covariates to correct for
#' @returns the table of the fit for the coefficient of interest
#' fit_ct <- do_differential_proportion_analysis(lpmcv2, cell_type_column = 'cell_type_safe', covariates=c('sex', 'age', 'day'))
do_differential_proportion_analysis <- function(seurat_object=NULL, metadata=NULL, cell_type_column='cell_type_final', variable_of_interest='timepoint.final', participant_column='assignment.final', covariates=c('lane', 'sex', 'age')) {
  # initialize variable
  meta.data <- NULL
  # extract metadata from object if required
  if (!is.null(metadata) & is.null(seurat_object)) {
    meta.data <- metadata
  }
  else if (is.null(metadata) & is.null(seurat_object)) {
    stop('metadata and Seurat object are both NULL, either needs to be non-NULL')
  }
  else if (!is.null(metadata) & !is.null(seurat_object)) {
    warning('both metadata and Seurat object are non-NULL, will use metadata')
    meta.data <- metadata
  }
  else if (is.null(metadata) & !is.null(seurat_object)) {
    meta.data <- seurat_object@meta.data
  }
  # collect the variables we need to combine
  combined_variables <- c(variable_of_interest, participant_column, covariates)
  # add as a column
  for (variable in combined_variables) {
    # as a new one if this is the first variable
    if (! ('combined_variable' %in% colnames(meta.data)) ) {
      meta.data[['combined_variable']] <- meta.data[[variable]]
    }
    # or add if we were already building
    else{
      meta.data[['combined_variable']] <- paste(meta.data[['combined_variable']], meta.data[[variable]], sep = '_')
    }
  }
  # convert to proportions
  props <- getTransformedProps(clusters = meta.data[[cell_type_column]], sample = meta.data[['combined_variable']])
  # get the unique combinations of covariates, variable of interest, and participant column
  pair_group <- unique(meta.data[, c('combined_variable', combined_variables)])
  # sort the pair group by this combination
  pair_group <- pair_group[order(pair_group[['combined_variable']]), ]
  # extract the 'pair', which is the donor
  pair <- pair_group[[participant_column]]
  # extract the 'group', which is the condition to compare
  group <- pair_group[[variable_of_interest]]
  # create the formula
  formula_to_use <- paste('~', variable_of_interest, '+', paste(covariates, collapse='+'), sep = '')
  des.tech <- model.matrix(as.formula(formula_to_use), data=pair_group)
  # create the correlations, with replicates (the donors)
  dupcor <- duplicateCorrelation(props$TransformedProps, design=des.tech,
                                 block=pair)
  # fit a regression model using limma
  fit1 <- lmFit(props$TransformedProps, design=des.tech, block=pair, 
                correlation=dupcor$consensus)
  # calculate statistics on the model
  fit1 <- eBayes(fit1)
  # get statistically significant proportions
  fit1_significant <- decideTests(fit1)
  # show summary
  summary(fit1_significant)
  # and specifically which ones
  topTable(fit1,coef=2, number = length(unique(meta.data[[cell_type_column]])))
  return(fit1)
}


####################
# Main Code        #
####################

# location of the metadata
stemi_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/cardio.integrated.20210301.metadata.tsv'
# where to place results
result_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/speckle/results/'

# extract the metadata
stemi_metadata <- read.table(stemi_metadata_loc, sep = '\t', row.names = 1, header = T)

# perform differential proportion analysis using the high resolution cell types
fit_ct <- do_differential_proportion_analysis(metadata = stemi_metadata, cell_type_column = 'cell_type_lowerres', covariates=c('chem', 'gender', 'age', 'assignment.final'))

# get statistically significant proportions
fit_ct_significant <- decideTests(fit_ct)
# show summary
summary(fit_ct_significant)
# and specifically which ones
topTable(fit_ct,coef=2, number = length(unique(stemi_metadata[['cell_type_lowerres']])))

# save result
fit_ct_df_t24h_t0 <- data.frame(topTable(fit_ct, coef=2, number = length(unique(stemi_metadata[['cell_type_lowerres']]))))
fit_ct_df_t24h_t0 <- cbind(data.frame(celltype=rownames(fit_ct_df_t24h_t0)), fit_ct_df_t24h_t0)
write.table(fit_ct_df_t24h_t0,
            paste(result_loc, 'timepoint_celltypelowerres_t24hBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# save result
fit_ct_df_t8w_t0 <- data.frame(topTable(fit_ct, coef=3, number = length(unique(stemi_metadata[['cell_type_lowerres']]))))
fit_ct_df_t8w_t0 <- cbind(data.frame(celltype=rownames(fit_ct_df_t8w_t0)), fit_ct_df_t8w_t0)
write.table(fit_ct_df_t8w_t0,
            paste(result_loc, 'timepoint_celltypelowerres_t8wBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# save result
fit_ct_df_C_t0 <- data.frame(topTable(fit_ct, coef=4, number = length(unique(stemi_metadata[['cell_type_lowerres']]))))
fit_ct_df_C_t0 <- cbind(data.frame(celltype=rownames(fit_ct_df_C_t0)), fit_ct_df_C_t0)
write.table(fit_ct_df_C_t0,
            paste(result_loc, 'timepoint_celltypelowerres_CBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# repeat for higher resolution cell types
fit_ct_hr <- do_differential_proportion_analysis(metadata = stemi_metadata, cell_type_column = 'cell_type', covariates=c('chem', 'gender', 'age'))
fit_ct_hr_significant <- decideTests(fit_ct_hr)
summary(fit_ct_hr_significant)
topTable(fit_ct_hr,coef=2, number = length(unique(stemi_metadata[['cell_type']])))
fit_ct_hr_df_t24h_t0 <- data.frame(topTable(fit_ct_hr, coef=2, number = length(unique(stemi_metadata[['cell_type']]))))
fit_ct_hr_df_t24h_t0 <- cbind(data.frame(celltype=rownames(fit_ct_hr_df_t24h_t0)), fit_ct_hr_df_t24h_t0)
fit_ct_hr_df_t8w_t0 <- data.frame(topTable(fit_ct_hr, coef=3, number = length(unique(stemi_metadata[['cell_type']]))))
fit_ct_hr_df_t8w_t0 <- cbind(data.frame(celltype=rownames(fit_ct_hr_df_t8w_t0)), fit_ct_hr_df_t8w_t0)
fit_ct_hr_df_C_t0 <- data.frame(topTable(fit_ct_hr, coef=4, number = length(unique(stemi_metadata[['cell_type']]))))
fit_ct_hr_df_C_t0 <- cbind(data.frame(celltype=rownames(fit_ct_hr_df_C_t0)), fit_ct_hr_df_C_t0)
write.table(fit_ct_hr_df_t24h_t0,
            paste(result_loc, 'timepoint_celltypehighres_t24hBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)
write.table(fit_ct_hr_df_t8w_t0,
            paste(result_loc, 'timepoint_celltypehighres_t8wBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)
write.table(fit_ct_hr_df_C_t0,
            paste(result_loc, 'timepoint_celltypehighres_CBaseline_20210301.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

