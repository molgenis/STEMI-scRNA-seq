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


####################
# main code        #
####################

# we need some more memory to do 
options(future.globals.maxSize = 120 * 1000 * 1024^2)

# the location of the objects
object_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
# specific objects
stemi_v2_loc <- paste(object_loc, 'stemi_final_wdemuxcorrectedassignments.rds', sep = '')
stemi_v3_loc <- paste(object_loc, 'v3_stemi_sct_idfixed_pca_SCTandRNAnorm.rds', sep = '')
hc_v2_loc <- paste(object_loc, 'HC_v2_20200616.rds', sep = '')
hc_v3_loc <- paste(object_loc, 'HC_v3_20200616.rds', sep = '')

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
stemi_v2@meta.data$assignment.final <- stemi_v2@meta.data$SNG.1ST
stemi_v2@meta.data$assignment_final <- NULL
stemi_v2@meta.data$timepoint.final <- stemi_v2@meta.data$timepoint.demux
stemi_v2@meta.data$assignment.ll <- stemi_v2@meta.data$assignment_ll
stemi_v2@meta.data$assignment_ll <- NULL
stemi_v2@meta.data$bare_barcode_lane <- rownames(stemi_v2@meta.data)
stemi_v2@meta.data$orig.ident <- 'stemi_v2'
stemi_v3@meta.data$orig.ident <- 'stemi_v3'
stemi_v3@meta.data$bare_barcode_lane <- rownames(stemi_v3@meta.data)


# add objects to list
cardio.list <- list(stemi_v2, stemi_v3, hc_v2, hc_v3)
# remove individual objects to clear memory
rm(stemi_v2)
rm(stemi_v3)
rm(hc_v2)
rm(hc_v3)
# select the features to use for integrating
cardio.features <- SelectIntegrationFeatures(object.list = cardio.list, nfeatures = 3000)
# perpare for the integration
cardio.list <- PrepSCTIntegration(object.list = cardio.list, anchor.features = cardio.features, verbose = T)
# get the anchors for integration
cardio.anchors <- FindIntegrationAnchors(object.list = cardio.list, normalization.method = "SCT", anchor.features = cardio.features, verbose = T)
# finally create the integrated object
cardio.integrated <- IntegrateData(anchorset = cardio.anchors, normalization.method = "SCT", verbose = T)
# clear memory by removing the list
rm(cardio.list)

# find variable features
cardio.integrated <- FindVariableFeatures(cardio.integrated)
# run PCA on the integrated object
cardio.integrated <- RunPCA(cardio.integrated)
# find the neighbours
cardio.integrated <- FindNeighbors(cardio.integrated, dims = 1:20)
# construct clusters from the neighbours
cardio.integrated <- FindClusters(cardio.integrated)
# perform UMAP dimension reduction for visualisation
cardio.integrated <- RunUMAP(cardio.integrated, dims = 1:20)

# grab our previously defined cell types
cell_types <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/metadata/cardio.integrated_cell_types_20200630.tsv', sep = '\t', header = T)
cardio.integrated <- AddMetaData(cardio.integrated, cell_types['cell_type'])
cardio.integrated <- AddMetaData(cardio.integrated, cell_types['cell_type_lowerres'])

# some may be undefined, we'll impute these, but we do need to keep track
cardio.integrated@meta.data$ct_was_imputed <- T
cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type), ]$ct_was_imputed <- F
# add imputed cell types
cardio.integrated <- add_imputed_meta_data(cardio.integrated, 'seurat_clusters', 'cell_type_lowerres', 'imputed_ct')
# for the missing cells, set this imputed cell type
cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$cell_type_lowerres), ]$cell_type_lowerres <- cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$cell_type_lowerres), ]$imputed_ct
# fix this little thing
cardio.integrated@meta.data[cardio.integrated@meta.data$orig.ident == 'stemi_v2', ]$chem <- 'V2'
cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$batch), ]$batch <- cardio.integrated@meta.data[is.na(cardio.integrated@meta.data$batch), ]$lane
saveRDS(cardio.integrated, paste(object_loc, 'cardio.integrated_20200625.rds', sep = ''))

# save our efforts
saveRDS(cardio.integrated, paste(object_loc, 'cardio.integrated_20200625.rds', sep = ''))
