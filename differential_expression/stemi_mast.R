####################
# libraries        #
####################

library(MAST)
library(Seurat)

####################
# Functions        #
####################

perform_mast <- function(seurat_object, output_loc, condition.1, condition.2, split.column = 'timepoint', assay = 'RNA', min.pct = 0.1){
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
    result <- FindMarkers(object = seurat_object, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, assay = assay)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t', col.names = 1)
  }
}

perform_mast_per_celltype <- function(seurat_object, output_loc, split.column = 'timepoint', cell.type.column = 'cell_type', assay = 'RNA', min.pct = 0.1){
  # go through the cell types
  for(cell_type in unique(seurat_object@meta.data[[cell.type.column]])){
    # make a more specific path
    output_loc_cell_type <- paste(output_loc, cell_type, sep = '')
    # grab the subset
    seurat_object_cell_type <- seurat_object[,seurat_object@meta.data[cell.type.column] == cell_type]
    # do variable feature detection, why Hilde?
    seurat_object <- FindVariableFeatures(object = seurat_object, mean.function = ExpMean, dispersion.function = LogVMR, nfeatures = 2000)
    # do the MAST
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'Baseline', condition.2 = 't24h', split.column = split.column, assay = assay, min.pct = min.pct)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 'Baseline', condition.2 = 't8w' ,split.column = split.column, assay = assay, min.pct = min.pct)
    perform_mast(seurat_object_cell_type, output_loc_cell_type, condition.1 = 't24h', condition.2 = 't8w' ,split.column = split.column, assay = assay, min.pct = min.pct)
  }
  # finally do a bulk analysis as well
  output_loc_bulk <- paste(output_loc, 'bulk', sep = '')
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'Baseline', condition.2 = 't24h', split.column = split.column, assay = assay, min.pct = min.pct)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 'Baseline', condition.2 = 't8w', split.column = split.column, assay = assay, min.pct = min.pct)
  perform_mast(seurat_object, output_loc_bulk, condition.1 = 't24h', condition.2 = 't8w', split.column = split.column, assay = assay, min.pct = min.pct)
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
#mast_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
mast_output_loc <- '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
# for a MAST comparison, also do only paired comparisons
mast_output_paired_loc_v2 <- paste(mast_output_loc, 'stemi_v2_paired/', sep = '')
mast_output_paired_loc_v3 <- paste(mast_output_loc, 'stemi_v3_paired/', sep = '')
mast_output_paired_lores_loc_v2 <- paste(mast_output_loc, 'stemi_v2_paired_lores/', sep = '')
mast_output_paired_lores_loc_v3 <- paste(mast_output_loc, 'stemi_v3_paired_lores/', sep = '')
mast_output_paired_loc_v2_rna <- paste(mast_output_paired_loc_v2, 'rna/', sep = '')
mast_output_paired_loc_v3_rna <- paste(mast_output_paired_loc_v3, 'rna/', sep = '')
mast_output_paired_loc_v2_sct <- paste(mast_output_paired_loc_v2, 'sct/', sep = '')
mast_output_paired_loc_v3_sct <- paste(mast_output_paired_loc_v3, 'sct/', sep = '')
mast_output_paired_lores_loc_v2_rna <- paste(mast_output_paired_lores_loc_v2, 'rna/', sep = '')
mast_output_paired_lores_loc_v3_rna <- paste(mast_output_paired_lores_loc_v3, 'rna/', sep = '')
mast_output_paired_lores_loc_v2_sct <- paste(mast_output_paired_lores_loc_v2, 'sct/', sep = '')
mast_output_paired_lores_loc_v3_sct <- paste(mast_output_paired_lores_loc_v3, 'sct/', sep = '')


# put in the work for v2
v2 <- readRDS(object_loc_v2)
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_sct, cell.type.column = 'cell_types', assay = 'SCT')
DefaultAssay(v2) <- 'RNA'
v2 <- NormalizeData(v2)
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_loc_v2_rna, cell.type.column = 'cell_types', assay = 'RNA')

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

# do the work again for with new cell types
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_sct, cell.type.column = 'cell_type_lowerres', assay = 'SCT')
perform_mast_per_celltype(seurat_object = v2, output_loc = mast_output_paired_lores_loc_v2_rna, cell.type.column = 'cell_type_lowerres', assay = 'RNA')

# clear up memory
rm(v2)
