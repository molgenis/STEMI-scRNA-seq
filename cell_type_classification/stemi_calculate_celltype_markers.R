#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_calculate_celltype_markers.R
# Function: calculate the celltype markers for the supplements
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(xlsx)


####################
# Functions        #
####################

list_of_dfs_to_excel <- function(list_of_dfs, excel_output_loc, remove_columns=NULL) {
  # create a new workbook
  wb = createWorkbook()
  # check each cell type
  for (item in names(list_of_dfs)) {
    # create the name for the sheet
    sheet_name <- item
    # create the sheet
    sheet = createSheet(wb, sheet_name)
    # fetch the dataframe
    item_result <- list_of_dfs[[item]]
    # again, replace labels if requested
    if (!is.null(remove_columns)) {
      item_result <- item_result[, setdiff(colnames(item_result), remove_columns)]
    }
    # add dataframe to sheet
    addDataFrame(item_result, sheet = sheet, startColumn=1, row.names=FALSE)
  }
  # write the result
  saveWorkbook(wb, excel_output_loc)
}

####################
# settings         #
####################

options(java.parameters = "-Xmx8000m")


####################
# debug code      #
####################

####################
# Main Code        #
####################

# location of the Seurat object
cardio_integrated_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'
# read the object
cardio_integrated <- readRDS(cardio_integrated_loc)

# calculate the marker genes for each cluster
Idents(cardio_integrated) <- 'seurat_clusters'
cardio_integrated_cluster_markers <- FindAllMarkers(cardio_integrated, test.use = 'wilcox')
cardio_integrated_cluster_markers <- cardio_integrated_cluster_markers[cardio_integrated_cluster_markers$p_val_adj < 0.05, ]

# calculate the marker genes for each lower resolution cell type
Idents(cardio_integrated) <- 'cell_type_lowerres'
cardio_integrated_cluster_low_celltypes <- FindAllMarkers(cardio_integrated, test.use = 'wilcox')
cardio_integrated_cluster_low_celltypes <- cardio_integrated_cluster_low_celltypes[cardio_integrated_cluster_low_celltypes$p_val_adj < 0.05, ]

# calculate the marker genes for each higher resolution cell type
Idents(cardio_integrated) <- 'cell_type'
cardio_integrated_cluster_high_celltypes <- FindAllMarkers(cardio_integrated, test.use = 'wilcox')
cardio_integrated_cluster_high_celltypes <- cardio_integrated_cluster_high_celltypes[cardio_integrated_cluster_high_celltypes$p_val_adj < 0.05, ]

# write results to tsv
write.table(cardio_integrated_cluster_markers, '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_cluster_marker_genes.tsv' , sep = '\t', row.names = F, col.names = T)
write.table(cardio_integrated_cluster_low_celltypes, '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_celltypelowerres_marker_genes.tsv' , sep = '\t', row.names = F, col.names = T)
write.table(cardio_integrated_cluster_high_celltypes, '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_celltypes_marker_genes.tsv' , sep = '\t', row.names = F, col.names = T)

# read all
cardio_integrated_cluster_markers <- read.table('/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_cluster_marker_genes.tsv' , sep = '\t', header = T)
cardio_integrated_cluster_low_celltypes <- read.table('/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_celltypelowerres_marker_genes.tsv' , sep = '\t',header = T)
cardio_integrated_cluster_high_celltypes <- read.table('/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_celltypes_marker_genes.tsv' , sep = '\t', header = T)
# put these markers in a list
markers_per_category <- list('clusters' = cardio_integrated_cluster_markers, 'cell_type' = cardio_integrated_cluster_high_celltypes, 'cell_type_lowerres' = cardio_integrated_cluster_low_celltypes)
# find a place to save these
marker_genes_save_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/cell_type_classification/stemi_marker_genes.xlsx'
# write to the file
list_of_dfs_to_excel(markers_per_category, marker_genes_save_loc)
