#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_compare_clustering_pcs.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)
library(ggplot2)
library(cowplot)


####################
# Functions        #
####################


create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL, legendless=F){
  # init the table
  confusion_table <- NULL
  # check each truth
  for(truth in unique(assignment_table[[truth_column]])){
    # get these truths
    truth_rows <- assignment_table[assignment_table[[truth_column]] == truth, ]
    # check now many have this truth
    this_truth_number <- nrow(truth_rows)
    # check what was predicted for these truths
    #for(prediction in unique(truth_rows[[prediction_column]])){
    for(prediction in unique(assignment_table[[prediction_column]])){
      # check the number of this prediction
      this_prediction_number <- nrow(truth_rows[truth_rows[[prediction_column]] == prediction, ])
      # init variable
      fraction <- NULL
      # we can only calculate a fraction if the result is not zero
      if(this_prediction_number > 0){
        # calculate the fraction
        fraction <- this_prediction_number / this_truth_number
      }
      # otherwise we just set it to zero
      else{
        fraction <- 0
      }
      # turn into row
      this_row <- data.frame(truth=c(truth), prediction=c(prediction), freq=c(fraction), stringsAsFactors = F)
      # add this entry to the dataframe
      if(is.null(confusion_table)){
        confusion_table <- this_row
      }
      else{
        confusion_table <- rbind(confusion_table, this_row)
      }
    }
  }
  # round the frequency off to a sensible cutoff
  confusion_table$freq <- round(confusion_table$freq, digits=2)
  # turn into plot
  p <- ggplot(data=confusion_table, aes(x=truth, y=prediction, fill=freq)) + geom_tile() + scale_fill_gradient(low='red', high='blue') + geom_text(aes(label=freq))
  # some options
  if(!is.null(truth_column_label)){
    p <- p + xlab(truth_column_label)
  }
  if(!is.null(prediction_column_label)){
    p <- p + ylab(prediction_column_label)
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}


# add metadata that is based on existing incomplete metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[[column_to_create]] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[!is.na(seurat_object@meta.data[[column_to_transform]]), column_to_transform])){
    # subset to get only this group
    seurat_group <- seurat_object[, !is.na(seurat_object@meta.data[[column_to_transform]]) & seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[!(is.na(seurat_group@meta.data[[column_to_reference]])) & seurat_group@meta.data[[column_to_reference]] == reference,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches", reference, correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[!is.na(seurat_object@meta.data[[column_to_transform]]) & seurat_object@meta.data[[column_to_transform]] == group, column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
  return(seurat_object)
}



####################
# Main Code        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

# location of the objects
cardio_integrated_loc <- '/groups/umcg-franke-scrna/tmp03/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'

# read the object
cardio_integrated <- readRDS(cardio_integrated_loc)
# update the object
cardio_integrated <- UpdateSeuratObject(cardio_integrated)
# set the correct assay
DefaultAssay(cardio_integrated) <- 'integrated'

# rerun the clustering an UMAP with 30 components. We set a seed this time for reproduciblity
cardio_integrated <- RunUMAP(cardio_integrated, dims = 1:30, return.model = T)
cardio_integrated <- FindNeighbors(cardio_integrated, dims = 1:30)
cardio_integrated <- FindClusters(cardio_integrated, resolution = 1.5)


# now we'll do the same with 50 PCs
cardio_integrated_pc50 <- RunUMAP(cardio_integrated, dims = 1:50, return.model = T)
cardio_integrated_pc50 <- FindNeighbors(cardio_integrated_pc50, dims = 1:50)
cardio_integrated_pc50 <- FindClusters(cardio_integrated_pc50, resolution = 1.5)

# first we'll compare the UMAPs
plot_grid(
  DimPlot(cardio_integrated, group.by = 'seurat_clusters'),
  DimPlot(cardio_integrated_pc50, group.by = 'seurat_clusters'),
  nrow = 1,
  ncol = 2
)

# now we'll check cluster membership
create_confusion_matrix(
  data.frame(
    pc30=cardio_integrated@meta.data[['seurat_clusters']], 
    pc50=cardio_integrated_pc50@meta.data[['seurat_clusters']]), 
  'pc30', 'pc50', 'leiden clusters 30 PCs', 'leiden clusters 50 PCs')

# do the same thing we did before
cardio_integrated@meta.data$cell_type <- cardio_integrated@meta.data$predicted.celltype.l2
cardio_integrated@meta.data$cell_type <- as.character(cardio_integrated@meta.data$cell_type)
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
cardio_integrated@meta.data[cardio_integrated@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
cardio_integrated@meta.data$cell_type <- gsub(' ', '_', cardio_integrated@meta.data$cell_type)
cardio_integrated <- add_imputed_meta_data(cardio_integrated, 'seurat_clusters', 'cell_type', 'ct_clus_prop')
cd4t <- c('Treg', 'CD4_Naive', 'CD4_TCM', 'CD4_TEM', 'CD4_CTL', 'CD4_Proliferating')
cd8t <- c('MAIT', 'CD8_Naive', 'CD8_TCM', 'CD8_TEM', 'CD8_Proliferating')
t_other <- c('dnT', 'gdT', 'ILC')
nk <- c('NKdim', 'NKbright', 'NK_Proliferating')
monocyte <- c('cMono', 'ncMono')
dc <- c('cDC1', 'cDC2', 'pDC', 'ASDC')
b <- c('B_naive', 'B_intermediate', 'B_memory')
cardio_integrated@meta.data$cell_type_lowerres <- cardio_integrated@meta.data$ct_clus_prop
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
#cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% t_other, ]$cell_type_lowerres <- 'T_other'
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% nk, ]$cell_type_lowerres <- 'NK'
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% dc, ]$cell_type_lowerres <- 'DC'
cardio_integrated@meta.data[cardio_integrated@meta.data$ct_clus_prop %in% b, ]$cell_type_lowerres <- 'B'

# do the same thing we did before
cardio_integrated_pc50@meta.data$cell_type <- cardio_integrated_pc50@meta.data$predicted.celltype.l2
cardio_integrated_pc50@meta.data$cell_type <- as.character(cardio_integrated_pc50@meta.data$cell_type)
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
cardio_integrated_pc50@meta.data$cell_type <- gsub(' ', '_', cardio_integrated_pc50@meta.data$cell_type)
cardio_integrated_pc50 <- add_imputed_meta_data(cardio_integrated_pc50, 'seurat_clusters', 'cell_type', 'ct_clus_prop')
cardio_integrated_pc50@meta.data$cell_type_lowerres <- cardio_integrated_pc50@meta.data$ct_clus_prop
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
#cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% t_other, ]$cell_type_lowerres <- 'T_other'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% nk, ]$cell_type_lowerres <- 'NK'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% dc, ]$cell_type_lowerres <- 'DC'
cardio_integrated_pc50@meta.data[cardio_integrated_pc50@meta.data$ct_clus_prop %in% b, ]$cell_type_lowerres <- 'B'

# now we'll check cluster membership
create_confusion_matrix(
  data.frame(
    pc30=cardio_integrated@meta.data[['cell_type_lowerres']], 
    pc50=cardio_integrated_pc50@meta.data[['cell_type_lowerres']]), 
  'pc30', 'pc50', 'cell types 30 PCs', 'cell types 50 PCs')
# saved as stemi_pc30_vs_pc50_celltypes