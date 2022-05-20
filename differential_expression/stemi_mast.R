####################
# libraries        #
####################

library(MAST)
library(Seurat)

####################
# Functions        #
####################

perform_mast <- function(seurat_object, output_loc, condition.1, condition.2, split.column = 'timepoint', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.25, latent.vars=NULL){
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[split.column], 'condition')
  Idents(seurat_object) <- seurat_object@meta.data$condition
  # create an output directory
  output_loc_final <- paste(output_loc, condition.1, condition.2, '.tsv', sep = '')
  # perform MAST if possible
  result <- NULL
  # we need to 'try' here, as MAST ungraciously exits if it cannot perform
  try({
    result <- FindMarkers(object = seurat_object, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, assay = assay, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t')
  }
}

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.25, latent.vars=NULL){
  DefaultAssay(seurat_object) <- assay
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do variable feature detection, why Hilde?
    seurat_object <- FindVariableFeatures(object = seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
    # do the MAST
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'Baseline', condition.2 = 't24h', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'Baseline', condition.2 = 't8w' ,split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 't24h', condition.2 = 't8w' ,split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = 'Baseline' ,split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = 't24h' ,split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'UT', condition.2 = 't8w' ,split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'Baseline', condition.2 = 't24h', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'Baseline', condition.2 = 't8w', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 't24h', condition.2 = 't8w', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = 'Baseline', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = 't24h', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'UT', condition.2 = 't8w', split.column = split.column, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
cardio_object_loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')

# DE output locations
mast_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
#mast_output_loc <- '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
# for a MAST comparison, also do only paired comparisons
mast_output_paired_loc_v2 <- paste(mast_output_loc, 'stemi_v2_paired_20210301/', sep = '')
mast_output_paired_loc_v3 <- paste(mast_output_loc, 'stemi_v3_paired_20210301/', sep = '')
mast_output_paired_lores_loc_v2 <- paste(mast_output_loc, 'stemi_v2_paired_lores_lfc01minpct01ncountrna_20210301/', sep = '')
mast_output_paired_lores_loc_v3 <- paste(mast_output_loc, 'stemi_v3_paired_lores_lfc01minpct01ncountrna_20210301/', sep = '')
mast_output_paired_loc_v2_rna <- paste(mast_output_paired_loc_v2, 'rna/', sep = '')
mast_output_paired_loc_v3_rna <- paste(mast_output_paired_loc_v3, 'rna/', sep = '')
mast_output_paired_loc_v2_sct <- paste(mast_output_paired_loc_v2, 'sct/', sep = '')
mast_output_paired_loc_v3_sct <- paste(mast_output_paired_loc_v3, 'sct/', sep = '')
mast_output_paired_lores_loc_v2_rna <- paste(mast_output_paired_lores_loc_v2, 'rna/', sep = '')
mast_output_paired_lores_loc_v3_rna <- paste(mast_output_paired_lores_loc_v3, 'rna/', sep = '')
mast_output_paired_lores_loc_v2_sct <- paste(mast_output_paired_lores_loc_v2, 'sct/', sep = '')
mast_output_paired_lores_loc_v3_sct <- paste(mast_output_paired_lores_loc_v3, 'sct/', sep = '')

# load object
cardio.integrated <- readRDS(cardio_object_loc)
# do the work for v2
cardio.chem2 <- subset(cardio.integrated, subset = chem == 'V2')
DefaultAssay(cardio.chem2) <- 'RNA'
cardio.chem2 <- NormalizeData(cardio.chem2)
perform_mast_per_celltype(cardio.chem2, mast_output_paired_lores_loc_v2_rna, split.column = 'timepoint.final', cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.1, latent.vars=c('nCount_RNA'))
rm(cardio.chem2)
# do the work for v3
cardio.chem3 <- subset(cardio.integrated, subset = chem == 'V3')
DefaultAssay(cardio.chem3) <- 'RNA'
cardio.chem3 <- NormalizeData(cardio.chem3)
perform_mast_per_celltype(cardio.chem3, mast_output_paired_lores_loc_v3_rna, split.column = 'timepoint.final', cell.type.column = 'cell_type_lowerres', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.1, latent.vars=c('nCount_RNA'))
rm(cardio.chem3)
