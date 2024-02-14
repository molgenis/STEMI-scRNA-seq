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
cardio_integrated_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'

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
