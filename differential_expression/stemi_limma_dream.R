
library('variancePartition')
library('edgeR')
library('BiocParallel')
library('Seurat')
library(Matrix)
library(Matrix.utils)

dreamer <- function(seurat_object, output_loc, aggregates=c('assignment.final', 'timepoint.final')){
  # grab the countmatrix
  countMatrix <- seurat_object@assays$RNA@counts
  # get the metadata
  metadata <- seurat_object@meta.data
  # create the groups to aggregate on, here it's on sample and timepoint
  groups <- metadata[, aggregates]
  # create an aggregated counts matrix
  aggregate_countMatrix <- t(aggregate.Matrix(t(countMatrix), groupings = groups, fun = 'sum'))
  # create aggregated metadata
  aggregate_metadata <- unique(metadata[, aggregates])
  # set the rownames of the aggregate metadata
  rownames_to_set_agg_metadata <- aggregate_metadata[[aggregates[1]]]
  for(i in 2:length(aggregates)){
    rownames_to_set_agg_metadata <- paste(rownames_to_set_agg_metadata, aggregate_metadata[[aggregates[i]]], sep='_')
  }
  rownames(aggregate_metadata) <- rownames_to_set_agg_metadata
  # set in the same order as the count matrix
  aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
  
  # filter genes by number of counts
  isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
  
  # Standard usage of limma/voom
  geneExpr = DGEList( aggregate_countMatrix[isexpr,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  # The variable to be tested must be a fixed effect
  form <- ~ timepoint.final + (1|assignment.final) 
  
  # estimate weights using linear mixed model of dream
  vobjDream = voomWithDreamWeights( geneExpr, form, aggregate_metadata )
  
  # Fit the dream model on each gene
  # By default, uses the Satterthwaite approximation for the hypothesis test
  #fitmm = dream( vobjDream, form, metadata )
  
  # use a contrast matrix
  # define and then cbind contrasts
  L1 = getContrast( vobjDream, form, aggregate_metadata, c("timepoint.finalt24h", "timepoint.finalt8w"))
  L2 = getContrast( vobjDream, form, aggregate_metadata, c("timepoint.finalt24h", "timepoint.finalBaseline"))
  L3 = getContrast( vobjDream, form, aggregate_metadata, c("timepoint.finalBaseline", "timepoint.finalt8w"))
  L = cbind(L1, L2, L3)     
  
  # fit both contrasts
  fit = dream( vobjDream, form, aggregate_metadata, L)
  
  # create list of variables
  vars <- list()
  vars[['form']] <- form
  vars[['vobjDream']] <- vobjDream
  vars[['L']] <- L
  vars[['fit']] <- fit
  
  saveRDS(vars, output_loc)
}

dream_pairwise <- function(seurat_object, output_loc, condition_combinations, aggregates=c('assignment.final', 'timepoint.final'), meta=F){
    # grab the countmatrix
    countMatrix <- seurat_object@assays$RNA@counts
    # get the metadata
    metadata <- seurat_object@meta.data
    # create the groups to aggregate on, here it's on sample and timepoint
    groups <- metadata[, aggregates]
    # create an aggregated counts matrix
    aggregate_countMatrix <- t(aggregate.Matrix(t(countMatrix), groupings = groups, fun = 'sum'))
    # create aggregated metadata
    aggregate_metadata <- unique(metadata[, aggregates])
    # set the rownames of the aggregate metadata
    rownames_to_set_agg_metadata <- aggregate_metadata[[aggregates[1]]]
    for(i in 2:length(aggregates)){
      rownames_to_set_agg_metadata <- paste(rownames_to_set_agg_metadata, aggregate_metadata[[aggregates[i]]], sep='_')
    }
    rownames(aggregate_metadata) <- rownames_to_set_agg_metadata
    # set in the same order as the count matrix
    aggregate_metadata <- aggregate_metadata[colnames(aggregate_countMatrix), ]
    # filter genes by number of counts
    isexpr = rowSums(cpm(aggregate_countMatrix)>0.1) >= 5
    
    # Standard usage of limma/voom
    geneExpr = DGEList( aggregate_countMatrix[isexpr,] )
    geneExpr = calcNormFactors( geneExpr )
    
    # Specify parallel processing parameters
    # this is used implicitly by dream() to run in parallel
    param = SnowParam(4, "SOCK", progressbar=TRUE)
    register(param)
    
    print(head(aggregate_metadata))
    
    # The variable to be tested must be a fixed effect
    form <- ~ 0 + timepoint.final + (1|assignment.final) 
    # do meta if requested
    if(meta){
      form <- ~ 0 + timepoint.final + chem + (1|assignment.final)
    }
    
    # estimate weights using linear mixed model of dream
    vobjDream = voomWithDreamWeights( geneExpr, form, aggregate_metadata )
    
    # do each combination
    for(combination_name in names(condition_combinations)){
      # grab the combination
      combination <- condition_combinations[[combination_name]]
      
      # define and then cbind contrasts
      L = getContrast( vobjDream, form, aggregate_metadata, c(paste('timepoint.final', combination[1], sep=''), paste('timepoint.final', combination[2], sep='')))
      
      # fit contrast
      fit = dream( vobjDream, form, aggregate_metadata, L)
      
      # grab the exact fit
      limma_result <- topTable(fit, coef=c('L1'), number=length(fit$F.p.value))
      
      # set an output location
      limma_output_loc <- paste(output_loc, combination_name, '.tsv', sep = '')
      # write the result
      write.table(limma_result, limma_output_loc, sep = '\t', row.names = T)
    }
    
}



# read command line arguments
args = commandArgs(trailingOnly=TRUE)
# get cell type and orig.ident
orig_ident <- args[1]
cell_type <- args[2]
# read the file
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20201209.rds')
# subset to the cell type and original identity
ct_and_oi <- cardio.integrated[, (cardio.integrated@meta.data$cell_type_lowerres == cell_type & cardio.integrated@meta.data$orig.ident == orig_ident)]
# clear some memory
rm(cardio.integrated)
# the output location
output_limma_base <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/limma/output/dream_contrast_'
output_limma <- paste(output_limma_base, orig_ident, '_', cell_type, '_pseudobulked.rds', sep='')
# call the method
dreamer(ct_and_oi, output_limma)


# make a list of condition combinations to test 1 vs 1
condition_combinations <- list()
ut_baseline <- c('UT', 'Baseline')
ut_t24h <- c('UT', 't24h')
ut_t8w <- c('UT', 't8w')
baseline_t24h <- c('Baseline', 't24h')
baseline_t8w <- c('Baseline', 't8w')
t24ht8w <- c('t24h', 't8w')
# add these to the list
condition_combinations[['ut_baseline']] <- ut_baseline
condition_combinations[['ut_t24h']] <- ut_t24h
condition_combinations[['ut_t8w']] <- ut_t8w
condition_combinations[['baseline_t24h']] <- baseline_t24h
condition_combinations[['baseline_t8w']] <- baseline_t8w
condition_combinations[['t24ht8w']] <- t24ht8w
# set the base path
limma_out_base_path <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/limma/output/dream_contrast_'
limma_out_ct_path <- paste(limma_out_base_path, cell_type, sep='')
# do the analysis
dream_pairwise(ct_and_oi, limma_out_ct_path, condition_combinations)
