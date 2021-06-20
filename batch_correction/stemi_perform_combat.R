# we need these libraries
library(Seurat)
library(sva)
library(Matrix)


runPhenograph <- function(seurat_object, k=30, dims=1:30){
  # get the assay
  assay_used <- seurat_object@reductions$pca@assay.used
  # calculate phenograph from the PCs
  seurat_object.pheno_out <- Rphenograph(as.matrix(seurat_object@reductions$pca@cell.embeddings[, dims]), k=k)
  # add to the metadata in the same way as the 'FindClusters' method
  seurat_object@meta.data$phenograph_clusters <- as.vector(unlist(membership(seurat_object.pheno_out.pheno_out[[2]])))
  seurat_object@meta.data[[paste(assay_used, 'pheno', k, sep = '_')]] <- as.vector(unlist(membership(seurat_object.pheno_out.pheno_out[[2]])))
  return(seurat_object)
}



# vector size needs to be increased
options(future.globals.maxSize = 370 * 1000 * 1024^2)
# load the file
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20210301.rds')
#cardio.integrated <- readRDS('/data/p287578/Cardiology/objects/cardio.integrated.20210301.rds')
# subset to STEMI data
cardio.stemi <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v2' | cardio.integrated@meta.data$orig.ident == 'stemi_v3']
# clear up memory
rm(cardio.integrated)
# set to SCT assay, as we'll correct those counts
DefaultAssay(cardio.stemi) <- 'SCT'
# perform combat
#corrected <- ComBat(dat=cardio.stemi@assays$SCT@counts, batch=as.factor(cardio.stemi@meta.data$chem))
#saveRDS('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/batch_correction/combat/cardio.stemi.20210301.combat_corrected.rds')
corrected_data <- ComBat(dat=cardio.stemi@assays$SCT@data, batch=as.factor(cardio.stemi@meta.data$chem))
#saveRDS(corrected_data, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/batch_correction/combat/cardio.stemi.20210301.combat_dataslot_corrected.rds')
saveRDS(corrected_data, '/data/p287578/Cardiology/objects/cardio.stemi.20210301.combat_dataslot_correction.rds')
# add combat data
cardio.stemi[['CBT']] <- CreateAssayObject(
  data=corrected_data[intersect(
    #rownames(cardio.stemi[, cardio.stemi@meta.data$orig.ident=='stemi_v2']@assays$SCT@counts), 
    #rownames(cardio.stemi[, cardio.stemi@meta.data$orig.ident=='stemi_v3']@assays$SCT@counts)
    rownames(cardio.stemi[, cardio.stemi@meta.data$orig.ident=='stemi_v2']@assays$SCT@data), 
    rownames(cardio.stemi[, cardio.stemi@meta.data$orig.ident=='stemi_v3']@assays$SCT@data)
  ), ]
)
# set the counts and scale.data as well
#cardio.stemi <- SetAssayData(cardio.stemi, 'CBT', 'counts', cardio.stemi@assays$SCT@counts[rownames( cardio.stemi@assays$CBT@data), ])
#cardio.stemi <- SetAssayData(cardio.stemi, 'CBT', 'scale.data', cardio.stemi@assays$CBT@data)
cardio.stemi$CBT@counts <-  cardio.stemi@assays$SCT@counts[rownames( cardio.stemi@assays$CBT@data), ]
cardio.stemi$CBT@scale.data <- cardio.stemi@assays$CBT@data
#
DefaultAssay(cardio.stemi) <- 'CBT'
# find variable features
cardio.stemi <- FindVariableFeatures(cardio.stemi)
# do PCA
cardio.stemi <- RunPCA(cardio.stemi, slot='data', assay='CBT')
# save the object
#saveRDS(cardio.stemi, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.stemi.20210611.combatcorrected.rds')
saveRDS(cardio.stemi, '/data/p287578/Cardiology/objects/cardio.stemi.20210611.combatdataslotcorrected.rds')
DefaultAssay(cardio.stemi) <- 'CBT'
cardio.stemi <- FindVariableFeatures(cardio.stemi)
cardio.stemi <- RunPCA(cardio.stemi)
cardio.stemi <- FindNeighbors(cardio.stemi, dims=1:30, assay='CBT', reduction='pca')
cardio.stemi <- FindClusters(cardio.stemi, resolution=1.2)
cardio.stemi <- RunUMAP(cardio.stemi, dims=1:30)
saveRDS(cardio.stemi, '/data/p287578/Cardiology/objects/cardio.stemi.20210611.combatdataslotcorrected_clustered.rds')
