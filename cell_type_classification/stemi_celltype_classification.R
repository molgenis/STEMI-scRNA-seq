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
# Main Code        #
####################

# locations of objects
object_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
object_version <- 'cardio.integrated.20201126_wazi.rds'
cardio.integrated_loc <- paste(object_loc, object_version, sep='')

# locations of the plots
dimplot_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/'
dimplot_ct_loc <- paste(dimplot_loc, 'cardio.integrated.20201209.cell_type.png', sep = '')
dimplot_ctl_loc <- paste(dimplot_loc, 'cardio.integrated.20201209.cell_type_lowerres.png', sep = '')

# read the object
cardio.integrated <- readRDS(cardio.integrated_loc)
# we will use the l2 predicted cell types (for the most part), so start by copying those
cardio.integrated@meta.data$cell_type <- cardio.integrated@meta.data$predicted.celltype.l2
# we'll convert it to a String, just to be sure that the next few steps are easier
cardio.integrated@meta.data$cell_type <- as.character(cardio.integrated@meta.data$cell_type)
# we'll rename some of the cell types so that they more closely follow our conventions
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
# spaces in variables is inconvenient in a lot of places, so we'll replace these with underscores
cardio.integrated@meta.data$cell_type <- gsub(' ', '_', cardio.integrated@meta.data$cell_type)

# start with the overwriting of some of the things we classified, there are likely no doublets, based on markers these should be ncMono
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'doublet', ]$cell_type <- 'ncMono'
# cluster 25 contains mostly erythrocytes, but based on expression, these are probably cMono, so overwriting those
cardio.integrated@meta.data[cardio.integrated@meta.data$seurat_clusters == 25, ]$cell_type <- 'cMono'
# for the rest of the erythrocytes, we will overwrite them with the greatest proportion cell type of the seurat cluster they are in, so first get that largest proportion
cardio.integrated <- add_imputed_meta_data(cardio.integrated, 'seurat_clusters', 'cell_type', 'ct_clus_prop')
# now overwrite for the erythrocytes
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'eryth', ]$cell_type <- cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'eryth', ]$ct_clus_prop

# we will define a lower resolution cell type as well, we need to create some groupings for this
cd4t <- c('Treg', 'CD4_Naive', 'CD4_TCM', 'CD4_TEM', 'CD4_CTL', 'CD4_Proliferating')
cd8t <- c('MAIT', 'CD8_Naive', 'CD8_TCM', 'CD8_TEM', 'CD8_Proliferating')
t_other <- c('dnT', 'gdT', 'ILC')
nk <- c('NKdim', 'NKbright', 'NK_Proliferating')
monocyte <- c('cMono', 'ncMono')
dc <- c('cDC1', 'cDC2', 'pDC')
b <- c('B_naive', 'B_intermediate', 'B_memory')
# add the new column by copying the higher res first, 
cardio.integrated@meta.data$cell_type_lowerres <- cardio.integrated@meta.data$cell_type
# in this new column, overwrite them to have the lower resolution
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% t_other, ]$cell_type_lowerres <- 'T_other'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% nk, ]$cell_type_lowerres <- 'NK'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% dc, ]$cell_type_lowerres <- 'DC'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type %in% b, ]$cell_type_lowerres <- 'B'

# change these into factors for convenience
cardio.integrated@meta.data$cell_type_lowerres <- as.factor(cardio.integrated@meta.data$cell_type_lowerres)
cardio.integrated@meta.data$cell_type <- as.factor(cardio.integrated@meta.data$cell_type)

# save this
saveRDS(cardio.integrated, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20201209.rds')

# see if this looks good
DimPlot(cardio.integrated, group.by = 'cell_type')
ggsave(dimplot_ct_loc, width=10, height=10)
DimPlot(cardio.integrated, group.by = 'cell_type_lowerres')
ggsave(dimplot_ctl_loc, width=10, height=10)
DimPlot(cardio.integrated, group.by = 'cell_type', label=T, label.size = 3, repel = T) + NoLegend()
ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.20201209.cell_type.nl.png', width=10, height=10)
DimPlot(cardio.integrated, group.by = 'cell_type_lowerres', label=T, label.size = 3, repel = T) + NoLegend()
ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.20201209.cell_type_lowerres.nl.png', width=10, height=10)

