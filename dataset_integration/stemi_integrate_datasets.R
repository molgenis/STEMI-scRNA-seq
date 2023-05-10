####################
# libraries        #
####################

library(Seurat)
library(ggplot2)

####################
# Functions        #
####################

# add metadata that is based on existing incomplete metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[column_to_create] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[[column_to_transform]])){
    # subset to get only this group
    seurat_group <- seurat_object[,seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[seurat_group@meta.data[[column_to_reference]] == reference & is.na(seurat_group@meta.data[[column_to_reference]]) == F,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[seurat_object@meta.data[[column_to_transform]] == group,][column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
  return(seurat_object)
}


add_gender_and_age <- function(seurat_object, gender_age_file_loc, assignment_column='assignment.final', column_to_add_age='age', column_to_add_gender='gender'){
  # grab the table from the file
  gender_age_mapping <- read.table(gender_age_file_loc, header = T, stringsAsFactors = F)
  # grab the barcodes and assignments
  assignments <- seurat_object@meta.data
  assignments$age <- NA
  assignments$gender <- NA
  # go through each row of the gender age file
  #apply(gender_age_mapping, 1, function(row){
  for(stemi_id in unique(gender_age_mapping$ID)){
    stemi_gender <- gender_age_mapping[gender_age_mapping$ID == stemi_id, 'gender']
    stemi_age <- gender_age_mapping[gender_age_mapping$ID == stemi_id, 'age']
    ll_id <- gender_age_mapping[gender_age_mapping$ID == stemi_id, 'll_match_id']
    ll_gender <- gender_age_mapping[gender_age_mapping$ID == stemi_id, 'll_match_gender']
    ll_age <- gender_age_mapping[gender_age_mapping$ID == stemi_id, 'll_match_age']
    # grab the values in the row
    #stemi_id <- row[['ID']]
    #stemi_gender <- row[['gender']]
    #stemi_age <- row[['age']]
    #ll_id <- row[['ll_match_id']]
    #ll_gender <- row[['ll_match_gender']]
    #ll_age <- row[['ll_match_age']]
    # add the age and gender
    if(nrow(assignments[assignments[[assignment_column]] == stemi_id, ]) > 0){
      assignments[assignments[[assignment_column]] == stemi_id, ]$age <- stemi_age
      assignments[assignments[[assignment_column]] == stemi_id, ]$gender <- stemi_gender
      
    }
    else{
      print(paste(stemi_id, 'not in object'))
    }
    if(nrow(assignments[assignments[[assignment_column]] == ll_id, ]) > 0){
      assignments[assignments[[assignment_column]] == ll_id, ]$age <- ll_age
      assignments[assignments[[assignment_column]] == ll_id, ]$gender <- ll_gender
    }
    else{
      print(paste(ll_id, 'not in object'))
    }
  }
  #})
  # add this info to the object
  seurat_object <- AddMetaData(seurat_object, assignments['age'], column_to_add_age)
  seurat_object <- AddMetaData(seurat_object, assignments['gender'], column_to_add_gender)
  return(seurat_object)
}

####################
# main code        #
####################

# we need some more memory to do 
options(future.globals.maxSize = 120 * 1000 * 1024^2)

# get age-gender file
age_gender_file_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/metadata/stemi_age_gender_match_wtestid.tsv'

# the location of the objects
object_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
# specific objects
stemi_v2_loc <- paste(object_loc, 'stemi_v2_normalized_samples_20201126.rds', sep = '')
stemi_v3_loc <- paste(object_loc, 'stemi_v3_normalized_samples_20201126.rds', sep = '')
hc_v2_loc <- paste(object_loc, 'HC_v2_normalized_samples_20201126.rds', sep = '')
hc_v3_loc <- paste(object_loc, 'HC_v3_normalized_samples_20201126.rds', sep = '')
cardio.integrated.loc <- paste(object_loc, 'cardio.integrated.20201126.rds', sep = '')

# read the objects
stemi_v2 <- readRDS(stemi_v2_loc)
stemi_v3 <- readRDS(stemi_v3_loc)
hc_v2 <- readRDS(hc_v2_loc)
hc_v3 <- readRDS(hc_v3_loc)

# make sure to set to the correct assay, also to make sure they have been SCTransformed
DefaultAssay(stemi_v2) <- 'SCT'
DefaultAssay(stemi_v3) <- 'SCT'
DefaultAssay(hc_v2) <- 'SCT'
DefaultAssay(hc_v3) <- 'SCT'

# fix some things in the objects
#stemi_v2@meta.data$assignment.final <- stemi_v2@meta.data$SNG.1ST
#stemi_v2@meta.data$assignment_final <- NULL
#stemi_v2@meta.data$timepoint.final <- stemi_v2@meta.data$timepoint.demux
#stemi_v2@meta.data$assignment.ll <- stemi_v2@meta.data$assignment_ll
#stemi_v2@meta.data$assignment_ll <- NULL
#stemi_v2@meta.data$bare_barcode_lane <- rownames(stemi_v2@meta.data)
stemi_v2@meta.data$orig.ident <- 'stemi_v2'
stemi_v3@meta.data$orig.ident <- 'stemi_v3'
#stemi_v3@meta.data$bare_barcode_lane <- rownames(stemi_v3@meta.data)

#v2 <- merge(hc_v2, stemi_v2)
#v2 <- SCTransform(v2, vars.to.regress = c('percent.mt'))
#v3 <- merge(hc_v3, stemi_v3)
#v3 <- SCTransform(v3, vars.to.regress = c('percent.mt'))
#cardio.list <- list(v2, v3)


# add objects to list
cardio.list <- list(stemi_v2, stemi_v3, hc_v2, hc_v3)
# remove individual objects to clear memory
rm(stemi_v2)
rm(stemi_v3)
rm(hc_v2)
rm(hc_v3)
# select the features to use for integrating
cardio.features <- SelectIntegrationFeatures(object.list = cardio.list, nfeatures = 2000)
# perpare for the integration
cardio.list <- PrepSCTIntegration(object.list = cardio.list, anchor.features = cardio.features, verbose = T)
# get the anchors for integration
cardio.anchors <- FindIntegrationAnchors(object.list = cardio.list, normalization.method = "SCT", anchor.features = cardio.features, verbose = T, k.filter = 200)
# finally create the integrated object
cardio.integrated <- IntegrateData(anchorset = cardio.anchors, normalization.method = "SCT", verbose = T)
# clear memory by removing the list
rm(cardio.list)

# find variable features
cardio.integrated <- FindVariableFeatures(cardio.integrated)
# run PCA on the integrated object
cardio.integrated <- RunPCA(cardio.integrated)
# find the neighbours
cardio.integrated <- FindNeighbors(cardio.integrated, dims = 1:30)
# construct clusters from the neighbours
cardio.integrated <- FindClusters(cardio.integrated, resolution = 1.5)
# perform UMAP dimension reduction for visualisation
cardio.integrated <- RunUMAP(cardio.integrated, dims = 1:30)
# make this a bit easier to plot
cardio.integrated@meta.data[cardio.integrated@meta.data$orig.ident == '1M_cells' & cardio.integrated@meta.data$chem == 'V2', ]$orig.ident <- 'hc_v2'
cardio.integrated@meta.data[cardio.integrated@meta.data$orig.ident == '1M_cells' & cardio.integrated@meta.data$chem == 'V3', ]$orig.ident <- 'hc_v3'
# do RNA norm
DefaultAssay(cardio.integrated) <- 'RNA'
cardio.integrated <- NormalizeData(cardio.integrated)
# add the age and gender
cardio.integrated <- add_gender_and_age(cardio.integrated, age_gender_file_loc)
# save the object
saveRDS(cardio.integrated.loc)

# grab our previously defined cell types
cell_types <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/metadata/cardio.integrated_cell_types_20200630.tsv', sep = '\t', header = T)
cardio.integrated <- AddMetaData(cardio.integrated, cell_types['cell_type'], 'cell_type.20200630')
cardio.integrated <- AddMetaData(cardio.integrated, cell_types['cell_type_lowerres'], 'cell_type_lowerres.20200630')

# read the azimuth imputed data
# azimuth_ct <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/seurat_multimodal/azimuth_cell_types_20201125.tsv', sep = '\t', row.names=1, header=T)
# cardio.integrated <- AddMetaData(cardio.integrated, azimuth_ct['predicted.celltype.l1.score'])
# cardio.integrated <- AddMetaData(cardio.integrated, azimuth_ct['predicted.celltype.l2.score'])
# cardio.integrated <- AddMetaData(cardio.integrated, azimuth_ct['predicted.celltype.l1'])
# cardio.integrated <- AddMetaData(cardio.integrated, azimuth_ct['predicted.celltype.l2'])
# FeaturePlot(cardio.integrated, features = c('predicted.celltype.l2.score'))

# # plot the imputed data
# DimPlot(cardio.integrated, label=T, label.size = 3, repel = T) + NoLegend()
# ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.30pcs.res1.2.clusters.20201125nl.png', width=10, height=10)
# DimPlot(cardio.integrated, group.by='cell_type.20200630', label=T, label.size = 3, repel = T) + NoLegend()
# ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.30pcs.res1.2.cell_type_20200630.20201125nl.png', width=10, height=10)
# DimPlot(cardio.integrated, group.by='cell_type_lowerres.20200630', label=T, label.size = 3, repel = T) + NoLegend()
# ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.30pcs.res1.2.cell_type_lowerres_20200630.20201125nl.png', width=10, height=10)
# DimPlot(cardio.integrated, group.by='predicted.celltype.l1', label=T, label.size = 3, repel = T) + NoLegend()
# ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.30pcs.res1.2.predicted.celltype.l1.20201125nl.png', width=10, height=10)
# DimPlot(cardio.integrated, group.by='predicted.celltype.l2', label=T, label.size = 3, repel = T) + NoLegend()
# ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.30pcs.res1.2.predicted.celltype.l2.20201125nl.png', width=10, height=10)
# # save the object
# saveRDS(cardio.integrated, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20201126_wazi.rds')

# # some may be undefined, we'll impute these, but we do need to keep track
# cardio.integrated@meta.data$ct_was_imputed <- T
# cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type), ]$ct_was_imputed <- F
# # add imputed cell types
# cardio.integrated <- add_imputed_meta_data(cardio.integrated, 'seurat_clusters', 'cell_type_lowerres', 'imputed_ct')
# # for the missing cells, set this imputed cell type
# cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$cell_type_lowerres), ]$cell_type_lowerres <- cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$cell_type_lowerres), ]$imputed_ct
# # fix this little thing
# #cardio.integrated@meta.data[cardio.integrated@meta.data$orig.ident == 'stemi_v2', ]$chem <- 'V2'
# cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$batch), ]$batch <- cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$batch), ]$lane


# # remove doublets mono 4
# cardio.integrated <- subset(cardio.integrated, subset = cell_type != 'mono 4')

# save our efforts
saveRDS(cardio.integrated, paste(object_loc, 'cardio.integrated_20200820.rds', sep = ''))
