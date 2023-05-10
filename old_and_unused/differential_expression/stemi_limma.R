####################
# libraries        #
####################

library(Seurat)
library(limma)

####################
# Functions        #
####################

perform_limma_trend <- function(seurat_object, output_loc, split.column = 'timepoint', assay = 'RNA', lfc = 1.2){
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[split.column], 'condition')
  # get expression from object
  expr <- as.matrix(GetAssayData(seurat_object))
  # get zero-expression genes
  bad <- which(rowSums(expr) == 0)
  # remove zero-expression genes
  expr <- expr[-bad,]
  # create model
  mm <- model.matrix(~0 + condition, data = seurat_object@meta.data)
  # create fit
  fit <- lmFit(expr, mm)
  # get DE with limma trend (lfc is significant logfold change)
  fit <- treat(fit, lfc=log2(lfc), trend = T)
  # get all the results
  result <- topTreat(fit, n = nrow(expr), sort.by = 'P')
  # save the result
  output_loc_final <- paste(output_loc, '.tsv', sep = '')
  write.table(result, output_loc_final, sep = '\t')
}

perform_limma_trend_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', lfc = 1.2){
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do the limma
    perform_limma_trend(seurat_object_cell_type, output_loc_cell_type, split.column = split.column, assay = assay, lfc = lfc)
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '_')
  perform_limma_trend(seurat_object, output_loc_bulk, split.column = split.column, assay = assay, lfc = lfc)
}

####################
# Main Code        #
####################

# object locations
object_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
#object_loc_v2 <- paste(object_loc, 'stemi_final.rds', sep = '')
object_loc_v2 <- paste(object_loc, 'stemi_final_wdemuxcorrectedassignments.rds', sep = '')
#object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200427.rds', sep = '')

# DE output locations
limma_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/limma/results/'
#limma_output_loc_v2 <- paste(limma_output_loc, 'stemi.final_v2/', sep = '')
limma_output_loc_v2 <- paste(limma_output_loc, 'stemi_v2/', sep = '')
limma_output_loc_v3 <- paste(limma_output_loc, 'stemi_v3/', sep = '')
limma_output_loc_v2_rna <- paste(limma_output_loc_v2, 'rna/', sep = '')
limma_output_loc_v3_rna <- paste(limma_output_loc_v3, 'rna/', sep = '')
limma_output_loc_v2_sct <- paste(limma_output_loc_v2, 'sct/', sep = '')
limma_output_loc_v3_sct <- paste(limma_output_loc_v3, 'sct/', sep = '')
# for a MAST comparison, also do only paired comparisons
limma_output_paired_loc_v2 <- paste(limma_output_loc, 'stemi_v2_paired/', sep = '')
limma_output_paired_loc_v3 <- paste(limma_output_loc, 'stemi_v3_paired/', sep = '')
limma_output_paired_loc_v2_rna <- paste(limma_output_paired_loc_v2, 'rna/', sep = '')
limma_output_paired_loc_v3_rna <- paste(limma_output_paired_loc_v3, 'rna/', sep = '')
limma_output_paired_loc_v2_sct <- paste(limma_output_paired_loc_v2, 'sct/', sep = '')
limma_output_paired_loc_v3_sct <- paste(limma_output_paired_loc_v3, 'sct/', sep = '')

#limma_output_loc_v2 <- paste(limma_output_loc, 'stemi.final_v2/', sep = '')
limma_output_lores_loc_v2 <- paste(limma_output_loc, 'stemi_v2_lores/', sep = '')
limma_output_lores_loc_v3 <- paste(limma_output_loc, 'stemi_v3_lores/', sep = '')
limma_output_lores_loc_v2_rna <- paste(limma_output_lores_loc_v2, 'rna/', sep = '')
limma_output_lores_loc_v3_rna <- paste(limma_output_lores_loc_v3, 'rna/', sep = '')
limma_output_lores_loc_v2_sct <- paste(limma_output_lores_loc_v2, 'sct/', sep = '')
limma_output_lores_loc_v3_sct <- paste(limma_output_lores_loc_v3, 'sct/', sep = '')
# for a MAST comparison, also do only paired comparisons
limma_output_paired_lores_loc_v2 <- paste(limma_output_loc, 'stemi_v2_paired_lores/', sep = '')
limma_output_paired_lores_loc_v3 <- paste(limma_output_loc, 'stemi_v3_paired_lores/', sep = '')
limma_output_paired_lores_loc_v2_rna <- paste(limma_output_paired_lores_loc_v2, 'rna/', sep = '')
limma_output_paired_lores_loc_v3_rna <- paste(limma_output_paired_lores_loc_v3, 'rna/', sep = '')
limma_output_paired_lores_loc_v2_sct <- paste(limma_output_paired_lores_loc_v2, 'sct/', sep = '')
limma_output_paired_lores_loc_v3_sct <- paste(limma_output_paired_lores_loc_v3, 'sct/', sep = '')

# put in the work for v2
v2 <- readRDS(object_loc_v2)
perform_limma_trend_per_celltype(seurat_object = v2, output_loc = limma_output_loc_v2_sct, cell.type.column = 'cell_types', assay = 'SCT')
DefaultAssay(v2) <- 'RNA'
v2 <- NormalizeData(v2)
perform_limma_trend_per_celltype(seurat_object = v2, output_loc = limma_output_loc_v2_rna, cell.type.column = 'cell_types', assay = 'RNA')
# now for the paired assignments
v2_baseline_t24h <- subset(v2, subset = timepoint == 'Baseline' | timepoint == 't24h')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t24h, output_loc = paste(limma_output_loc_v2_rna, 'baseline_t24h_', sep = ''), cell.type.column = 'cell_types', assay = 'RNA')
v2_baseline_t8w <- subset(v2, subset = timepoint == 'Baseline' | timepoint == 't8w')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t8w, output_loc = paste(limma_output_loc_v2_rna, 'baseline_t8w_', sep = ''), cell.type.column = 'cell_types', assay = 'RNA')
v2_t24h_t8w <- subset(v2, subset = timepoint == 't24h' | timepoint == 't8w')
perform_limma_trend_per_celltype(seurat_object = v2_t24h_t8w, output_loc = paste(limma_output_loc_v2_rna, 't24h_t8w_', sep = ''), cell.type.column = 'cell_types', assay = 'RNA')

perform_limma_trend_per_celltype(seurat_object = v2_baseline_t24h, output_loc = paste(limma_output_loc_v2_sct, 'baseline_t24h_', sep = ''), cell.type.column = 'cell_types', assay = 'SCT')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t8w, output_loc = paste(limma_output_loc_v2_sct, 'baseline_t8w_', sep = ''), cell.type.column = 'cell_types', assay = 'SCT')
perform_limma_trend_per_celltype(seurat_object = v2_t24h_t8w, output_loc = paste(limma_output_loc_v2_sct, 't24h_t8w_', sep = ''), cell.type.column = 'cell_types', assay = 'SCT')

# downsample cell types
v2@meta.data$cell_type_lowerres <- v2@meta.data$cell_types
levels(v2@meta.data$cell_type_lowerres) <- c(levels(v2@meta.data$cell_type_lowerres), 'NK', 'DC', 'monocyte', 'CD4T', 'CD8T')
v2@meta.data[v2@meta.data$cell_type_lowerres == 'NKdim', ]$cell_type_lowerres <- 'NK'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'NKbright', ]$cell_type_lowerres <- 'NK'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'cMonocytes', ]$cell_type_lowerres <- 'monocyte'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'ncMonocytes', ]$cell_type_lowerres <- 'monocyte'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'mDCs', ]$cell_type_lowerres <- 'DC'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'pDCs', ]$cell_type_lowerres <- 'DC'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'Naive_CD8T', ]$cell_type_lowerres <- 'CD8T'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'Memory_CD8T', ]$cell_type_lowerres <- 'CD8T'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'Naive_CD4T', ]$cell_type_lowerres <- 'CD4T'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'Memory_CD4T', ]$cell_type_lowerres <- 'CD4T'
v2@meta.data[v2@meta.data$cell_type_lowerres == 'PlasmaB', ]$cell_type_lowerres <- 'B'

# do the whole shebang again
perform_limma_trend_per_celltype(seurat_object = v2, output_loc = limma_output_lores_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_celltype(seurat_object = v2, output_loc = limma_output_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')

# now for the paired assignments
v2_baseline_t24h <- subset(v2, subset = timepoint == 'Baseline' | timepoint == 't24h')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t24h, output_loc = paste(limma_output_lores_loc_v2_rna, 'baseline_t24h_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'RNA')
v2_baseline_t8w <- subset(v2, subset = timepoint == 'Baseline' | timepoint == 't8w')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t8w, output_loc = paste(limma_output_lores_loc_v2_rna, 'baseline_t8w_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'RNA')
v2_t24h_t8w <- subset(v2, subset = timepoint == 't24h' | timepoint == 't8w')
perform_limma_trend_per_celltype(seurat_object = v2_t24h_t8w, output_loc = paste(limma_output_lores_loc_v2_rna, 't24h_t8w_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'RNA')

perform_limma_trend_per_celltype(seurat_object = v2_baseline_t24h, output_loc = paste(limma_output_lores_loc_v2_sct, 'baseline_t24h_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_celltype(seurat_object = v2_baseline_t8w, output_loc = paste(limma_output_lores_loc_v2_sct, 'baseline_t8w_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_limma_trend_per_celltype(seurat_object = v2_t24h_t8w, output_loc = paste(limma_output_lores_loc_v2_sct, 't24h_t8w_', sep = ''), cell.type.column = 'cell_type_lowerres', assay = 'SCT')

# clear up memory
rm(v2)
