library(Seurat)
library(SeuratDisk)
library(ggplot2)
library(patchwork)

create_confusion_matrix <- function(assignment_table, truth_column, prediction_column, truth_column_label=NULL, prediction_column_label=NULL){
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
  return(p)
}



# load the reference Seurat object
reference <- LoadH5Seurat('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/cell-type-classifying/seurat_multimodal/references/pbmc_multimodal.h5seurat')


# load the query object
stemi_v3 <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v3_normalized_samples_20201126.rds')
# set to the assay used in the reference
DefaultAssay(stemi_v3) <- 'SCT'

# find transfer anchors between the reference and the query, the query is your dataset
anchors.stemi_v3 <- FindTransferAnchors(
  reference = reference,
  query = stemi_v3,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
stemi_v3 <- MapQuery(
  anchorset = anchors.stemi_v3,
  query = stemi_v3,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# save the new object
saveRDS(stemi_v3, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v3_normalized_samples_20210512_wazi.rds')

# load the query object
stemi_v2 <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v2_normalized_samples_20201126.rds')
# set to the assay used in the reference
DefaultAssay(stemi_v2) <- 'SCT'

# find transfer anchors between the reference and the query, the query is your dataset
anchors.stemi_v2 <- FindTransferAnchors(
  reference = reference,
  query = stemi_v2,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
stemi_v2 <- MapQuery(
  anchorset = anchors.stemi_v2,
  query = stemi_v2,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# save the new object
saveRDS(stemi_v2, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v2_normalized_samples_20210512_wazi.rds')

# load the query object
hc_v3 <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/HC_v3_normalized_samples_20201126.rds')
# set to the assay used in the reference
DefaultAssay(hc_v3) <- 'SCT'

# find transfer anchors between the reference and the query, the query is your dataset
anchors.hc_v3 <- FindTransferAnchors(
  reference = reference,
  query = hc_v3,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
hc_v3 <- MapQuery(
  anchorset = anchors.hc_v3,
  query = hc_v3,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# save the new object
saveRDS(hc_v3, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/HC_v3_normalized_samples_20210512_wazi.rds')


# load the query object
hc_v2 <- readRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/HC_v2_normalized_samples_20201126.rds')
# set to the assay used in the reference
DefaultAssay(hc_v2) <- 'SCT'

# find transfer anchors between the reference and the query, the query is your dataset
anchors.hc_v2 <- FindTransferAnchors(
  reference = reference,
  query = hc_v2,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE
)

# map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
hc_v2 <- MapQuery(
  anchorset = anchors.hc_v2,
  query = hc_v2,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca",
  reduction.model = "wnn.umap"
)

# save the new object
saveRDS(hc_v2, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/HC_v2_normalized_samples_20210512_wazi.rds')

# summarize this new classification
predictions_hc_v2 <- hc_v2@meta.data[, c('predicted.celltype.l1', 'predicted.celltype.l1.score', 'predicted.celltype.l2', 'predicted.celltype.l2.score')]
predictions_hc_v3 <- hc_v3@meta.data[, c('predicted.celltype.l1', 'predicted.celltype.l1.score', 'predicted.celltype.l2', 'predicted.celltype.l2.score')]
predictions_stemi_v2 <- stemi_v2@meta.data[, c('predicted.celltype.l1', 'predicted.celltype.l1.score', 'predicted.celltype.l2', 'predicted.celltype.l2.score')]
predictions_stemi_v3 <- stemi_v3@meta.data[, c('predicted.celltype.l1', 'predicted.celltype.l1.score', 'predicted.celltype.l2', 'predicted.celltype.l2.score')]
# add into one table
predictions_all <- rbind(predictions_stemi_v2, predictions_stemi_v3, predictions_hc_v2, predictions_hc_v3)
# add the barcode as an explicit column
predictions_all$barcode <- rownames(predictions_all)
# change the order
predictions_all <- predictions_all[, c('barcode', 'predicted.celltype.l1', 'predicted.celltype.l1.score', 'predicted.celltype.l2', 'predicted.celltype.l2.score')]
# write the table
write.table(predictions_all, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/seurat_multimodal/azimuth_cell_types_20210521.tsv', sep = '\t', col.names = T, row.names = F)
# clear some memory
rm(hc_v2)
rm(hc_v3)
rm(stemi_v2)
rm(stemi_v3)

# # check concordance with previous classification
# cardio.integrated <- AddMetaData(cardio.integrated, predictions_all['predicted.celltype.l1'], 'predicted.celltype.l1.v4')
# create_confusion_matrix(cardio.integrated@meta.data[, c('predicted.celltype.l1', 'predicted.celltype.l1.v4')], 'predicted.celltype.l1', 'predicted.celltype.l1.v4', 'prerelease', 'final release')
# ggsave('cardio.integrated.azimuth_l1_confusion_matrix_pre_and_final_release.pdf', width=20, height=20)
# cardio.integrated <- AddMetaData(cardio.integrated, predictions_all['predicted.celltype.l2'], 'predicted.celltype.l2.v4')
# create_confusion_matrix(cardio.integrated@meta.data[, c('predicted.celltype.l2', 'predicted.celltype.l2.v4')], 'predicted.celltype.l2', 'predicted.celltype.l2.v4', 'prerelease', 'final release')
# ggsave('cardio.integrated.azimuth_l2_confusion_matrix_pre_and_final_release.pdf', width=20, height=20)


# try with the integrated object as well
storage_read_part <- 'tmp01'
storage_write_part <- 'tmp01'

# object locations
object_loc <- paste('/groups/umcg-wijmenga/', storage_read_part, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/', sep = '')
cardio_object_unclassificed_loc <- paste(object_loc, 'cardio.integrated_20200820.rds', sep = '')
# read the object
cardio_integrated <- readRDS(cardio_object_unclassificed_loc)
# add the celltype classification
cardio.integrated <- AddMetaData(cardio.integrated, predictions_all['predicted.celltype.l1'], 'predicted.celltype.l1')
cardio.integrated <- AddMetaData(cardio.integrated, predictions_all['predicted.celltype.l2'], 'predicted.celltype.l2')
# get path of new object
cardio_object_classificed_loc <- paste(object_loc, 'cardio.integrated.20201126_wazi.rds', sep = '')
# save result
saveRDS(cardio_integrated, cardio_object_classificed_loc)
