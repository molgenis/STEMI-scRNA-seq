
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
  # set the rownames to be the same as the countsMatrix aggregate colnames, should go correctly if data is valid
  rownames(aggregate_metadata) <- colnames(aggregate_countMatrix)
  
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
  L1 = getContrast( vobjDream, form, metadata_metadata, c("timepoint.finalt24h", "timepoint.finalt8w"))
  L2 = getContrast( vobjDream, form, metadata_metadata, c("timepoint.finalt24h", "timepoint.finalBaseline"))
  L3 = getContrast( vobjDream, form, metadata_metadata, c("timepoint.finalBaseline", "timepoint.finalt8w"))
  L = cbind(L1, L2, L3)     
  
  # fit both contrasts
  fit = dream( vobjDream, form, metadata_metadata, L)
  
  # create list of variables
  vars <- list()
  vars[['form']] <- form
  vars[['vobjDream']] <- vobjDream
  vars[['L']] <- L
  vars[['fit']] <- fit
  
  saveRDS(vars, output_loc)
}


# read command line arguments
args = commandArgs(trailingOnly=TRUE)
# get cell type and orig.ident
orig_ident <- args[1]
cell_type <- args[2]
# read the file
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated_20200820.rds')
# subset to the cell type and original identity
ct_and_oi <- cardio.integrated[, (cardio.integrated@meta.data$cell_type_lowerres == cell_type & cardio.integrated@meta.data$orig.ident == orig_ident)]
# clear some memory
rm(cardio.integrated)
# the output location
output_limma_base <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/limma/output/dream_contrast_'
output_limma <- paste(output_limma_base, orig_ident, '_', cell_type, '.rds')
# call the method
dreamer(ct_and_oi, output_limma)

pilot3 <- readRDS('/data/scRNA/Seurat/objects/pilot3_seurat3_200420_phenofixed_ctfixed.rds')