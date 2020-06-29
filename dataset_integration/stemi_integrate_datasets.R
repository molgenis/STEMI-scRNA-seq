####################
# libraries        #
####################

library(Seurat)
library(ggplot2)

####################
# Functions        #
####################



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

# save our efforts
saveRDS(cardio.integrated, paste(object_loc, 'cardio.integrated_20200625.rds', sep = ''))
