#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_mark_communicating_monocytes.R
# Function: mark the monocytes that are communicating
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

get_interactions_from_table <- function(interaction_table, cutoff=0.05){
  # we have a vector in which we will store the interactions
  interactions <- c()
  # check each of the downstream genes
  downstream_genes <- rownames(interaction_table)
  # check each of the ligands
  ligands <- colnames(interaction_table)
  # now check each combination
  for(downstream_gene in downstream_genes){
    for(ligand in ligands){
      # check the interaction
      interaction <- interaction_table[downstream_gene, ligand]
      # check if the interaction is large enough
      if(interaction > cutoff){
        # add it to the list of interactions
        interactions <- c(interactions, paste(ligand, downstream_gene, sep = '_'))
      }
    }
  }
  return(interactions)
}


combine_interactions_per_pathways <- function(interactions_per_pathway_1, interactions_per_pathway_2, cutoff=0.05, by_receptor=T){
  # store the interactions per pathway
  interactions_per_pathway_combined <- list()
  # check which pathways we can use
  common_pathways <- intersect(names(interactions_per_pathway_1), names(interactions_per_pathway_2))
  # check these
  for(pathway in common_pathways){
    # we will save the result
    interaction_numbers <- NULL
    # grab for that specific pathway
    nichenet_list_pathway_1 <- interactions_per_pathway_1[[pathway]]
    nichenet_list_pathway_2 <- interactions_per_pathway_2[[pathway]]
    # check which cell types we can do
    cell_types_senders_to_check <- intersect(names(nichenet_list_pathway_1), names(nichenet_list_pathway_2))
    # check each of these cell types
    for(cell_type_sending in cell_types_senders_to_check){
      # subset to that
      nichenet_receivers_list_1 <- nichenet_list_pathway_1[[cell_type_sending]]
      nichenet_receivers_list_2 <- nichenet_list_pathway_2[[cell_type_sending]]
      # now check which of these overlap
      cell_types_receivers_to_check <- intersect(names(nichenet_receivers_list_1), names(nichenet_receivers_list_2))
      # check each of these cell types
      for(cell_type_receiving in cell_types_receivers_to_check){
        # grab the interactions
        interactions_1 <- NULL
        interactions_2 <- NULL
        if(by_receptor){
          interactions_1 <- nichenet_receivers_list_1[[cell_type_receiving]][['ligand_receptor_network']]
          interactions_2 <- nichenet_receivers_list_2[[cell_type_receiving]][['ligand_receptor_network']]
        }
        else{
          interactions_1 <- nichenet_receivers_list_1[[cell_type_receiving]][['ligand_target']]
          interactions_2 <- nichenet_receivers_list_2[[cell_type_receiving]][['ligand_target']]
        }
        # get the interactions from the tables
        interactions_genes_1 <- get_interactions_from_table(interactions_1, cutoff = cutoff)
        interactions_genes_2 <- get_interactions_from_table(interactions_2, cutoff = cutoff)
        # get the outer join of the genes
        interactions_genes <- unique(c(interactions_genes_1, interactions_genes_2))
        # count the number of interactions
        interaction_number <- length(interactions_genes)
        # create a row
        interaction_row <- data.frame(from=c(cell_type_sending), to=c(cell_type_receiving), number=c(interaction_number))
        # if this is the first round through, we will initialize the dataframe
        if(is.null(interaction_numbers)){
          interaction_numbers <- interaction_row
        }
        else{
          interaction_numbers <- rbind(interaction_numbers, interaction_row)
        }
      }
    }
    interactions_per_pathway_combined[[pathway]] <- interaction_numbers
  }
  return(interactions_per_pathway_combined)
}


get_downstream_genes_from_celltype <- function(nichenet_receiver, cutoff=0.05){
  # collect all the downstream genes
  downstream_affected_genes <- c()
  # check each cell type
  for(sending_celltype in names(nichenet_receiver)){
    # get the output
    sending_cells <- nichenet_receiver[[sending_celltype]]
    # check if there is a ligand-target
    if('ligand_target' %in% names(sending_cells)){
      # get the target table
      target_table <- sending_cells[['ligand_target']]
      # check if there are any significant interactions
      nr_sig_per_gene <- apply(target_table, 2, function(x){
        sum(x > 0.05)
      })
      # get the genes where there was one or more significant interactions
      genes_downstream_sig <- names(nr_sig_per_gene)[nr_sig_per_gene > 0]
      # add to the list
      downstream_affected_genes <- c(downstream_affected_genes, genes_downstream_sig)
    }
  }
  # filter what might be there multiple times
  downstream_affected_genes <- unique(downstream_affected_genes)
  return(downstream_affected_genes)
}


####################
# Main Code        #
####################

# location of the objects
objects.loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
cardio.integrated.loc <- paste(objects_loc, 'cardio.integrated.20210301.rds', sep = '')
# read the object
cardio.integrated <- readRDS(cardio.integrated.loc)

# location of the nichenet output
nichenet.v2.Baseline.vs.t24h.loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v2_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted_perpathway_20220621.rds'
nichenet.v3.Baseline.vs.t24h.loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v3_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted_perpathway_20220621.rds'
nichenet.v2.Baseline.vs.t8w.loc<- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v2_Baseline_vs_t8w_nichenet_onlymajors_perct_omni_unweighted_perpathway_20220621.rds'
nichenet.v3.Baseline.vs.t8w.loc<- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v3_Baseline_vs_t8w_nichenet_onlymajors_perct_omni_unweighted_perpathway_20220621.rds'
# read the object
nichenet.v2.Baseline.vs.t24h <- readRDS(nichenet.v2.Baseline.vs.t24h.loc)
nichenet.v3.Baseline.vs.t24h <- readRDS(nichenet.v3.Baseline.vs.t24h.loc)
nichenet.v2.Baseline.vs.t8w <- readRDS(nichenet.v2.Baseline.vs.t8w.loc)
nichenet.v3.Baseline.vs.t8w <- readRDS(nichenet.v3.Baseline.vs.t8w.loc)

# get the il1 table for monocytes
il1.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`il1-pathway-reactome`$monocyte)
il1.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`il1-pathway-reactome`$monocyte)
# now IL4-IL13 signalling
il4.13.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`il4il13-pathway-reactome`$monocyte)
il4.13.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`il4il13-pathway-reactome`$monocyte)
il10.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`il10-pathway-reactome`$monocyte)
il10.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`il10-pathway-reactome`$monocyte)
il17.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`il17-pathway-reactome`$monocyte)
il17.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`il17-pathway-reactome`$monocyte)
il6.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`il6-pathway-reactome`$monocyte)
il6.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`il6-pathway-reactome`$monocyte)
il17.kegg.v2.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t24h$`KEGG IL-17 signaling pathway`$monocyte)
il17.kegg.v3.Baseline.vs.t24h.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t24h$`KEGG IL-17 signaling pathway`$monocyte)
# also for t8w
il1.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`il1-pathway-reactome`$monocyte)
il1.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`il1-pathway-reactome`$monocyte)
il4.13.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`il4il13-pathway-reactome`$monocyte)
il4.13.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`il4il13-pathway-reactome`$monocyte)
il10.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`il10-pathway-reactome`$monocyte)
il10.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`il10-pathway-reactome`$monocyte)
il17.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`il17-pathway-reactome`$monocyte)
il17.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`il17-pathway-reactome`$monocyte)
il6.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`il6-pathway-reactome`$monocyte)
il6.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`il6-pathway-reactome`$monocyte)
il17.kegg.v2.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v2.Baseline.vs.t8w$`KEGG IL-17 signaling pathway`$monocyte)
il17.kegg.v3.Baseline.vs.t8w.mono.targets <- get_downstream_genes_from_celltype(nichenet.v3.Baseline.vs.t8w$`KEGG IL-17 signaling pathway`$monocyte)


# extract STEMIv2 monocytes
stemi.v2.monocytes <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v2' & cardio.integrated@meta.data$cell_type_lowerres == 'monocyte']
# do the PCA/UMAP/etc
DefaultAssay(stemi.v2.monocytes) <- 'SCT'
stemi.v2.monocytes <- FindVariableFeatures(stemi.v2.monocytes)
# set the seed for reproducability
set.seed(7777)
stemi.v2.monocytes <- RunPCA(stemi.v2.monocytes)
stemi.v2.monocytes <- RunUMAP(stemi.v2.monocytes, dim = 1:30, reduction = 'pca')
stemi.v2.monocytes <- FindNeighbors(stemi.v2.monocytes)
stemi.v2.monocytes <- FindClusters(stemi.v2.monocytes)
# extract STEMIv3 monocytes
stemi.v3.monocytes <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v3' & cardio.integrated@meta.data$cell_type_lowerres == 'monocyte']
# do the PCA/UMAP/etc
DefaultAssay(stemi.v3.monocytes) <- 'SCT'
stemi.v3.monocytes <- FindVariableFeatures(stemi.v3.monocytes)
stemi.v3.monocytes <- RunPCA(stemi.v3.monocytes)
stemi.v3.monocytes <- RunUMAP(stemi.v3.monocytes, dim = 1:30, reduction = 'pca')
stemi.v3.monocytes <- FindNeighbors(stemi.v3.monocytes)
stemi.v3.monocytes <- FindClusters(stemi.v3.monocytes)
# check locations of timepoints
stemi.v2.monocytes.dimplot.timepoint <- DimPlot(stemi.v2.monocytes, group.by = 'timepoint.final')
stemi.v3.monocytes.dimplot.timepoint <- DimPlot(stemi.v3.monocytes, group.by = 'timepoint.final')
# check the expression of the downstream genes
stemi.v2.monocytes.dimplot.downstream.features.il4.13 <- FeaturePlot(stemi.v2.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets))
stemi.v3.monocytes.dimplot.downstream.features.il4.13 <- FeaturePlot(stemi.v3.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets))
# save the plots somewhere
features_plot_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/downstream_features_plots/'
ggsave(plot = stemi.v2.monocytes.dimplot.downstream.features.il4.13, filename = 'stemi_v2_monocytes_downstream_features_il4_il13.pdf', path = features_plot_loc, width = 5, height = 5)
ggsave(plot = stemi.v3.monocytes.dimplot.downstream.features.il4.13, filename = 'stemi_v2_monocytes_downstream_features_il4_il13.pdf', path = features_plot_loc, width = 5, height = 5)
stemi.v2.monocytes <- AddModuleScore(stemi.v2.monocytes, features = list('il4-13'=c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets)), name = 'il4-13')
stemi.v3.monocytes <- AddModuleScore(stemi.v3.monocytes, features = list('il4-13'=c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets)), name = 'il4-13')
FeaturePlot(stemi.v2.monocytes, features = c('il4.131'))
FeaturePlot(stemi.v3.monocytes, features = c('il4.131'))
# also t8w
stemi.v2.monocytes <- AddModuleScore(stemi.v2.monocytes, features = list('il4-13-t8w'=c(il4.13.v2.Baseline.vs.t8w.mono.targets, il4.13.v3.Baseline.vs.t8w.mono.targets)), name = 'il4-13-t8w')
stemi.v3.monocytes <- AddModuleScore(stemi.v3.monocytes, features = list('il4-13-t8w'=c(il4.13.v2.Baseline.vs.t8w.mono.targets, il4.13.v3.Baseline.vs.t8w.mono.targets)), name = 'il4-13-t8w')
FeaturePlot(stemi.v2.monocytes, features = c('il4.13.t8w1'))
FeaturePlot(stemi.v3.monocytes, features = c('il4.13.t8w1'))

# get the average expression per timepoint
stemi.v2.monocytes.exp <- AverageExpression(stemi.v2.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets), assays = 'SCT', group.by = 'timepoint.final')
stemi.v3.monocytes.exp <- AverageExpression(stemi.v3.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets), assays = 'SCT', group.by = 'timepoint.final')
stemi.v2.monocytes.clust.exp <- AverageExpression(stemi.v2.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets), assays = 'SCT', group.by = 'seurat_clusters')
stemi.v3.monocytes.clust.exp <- AverageExpression(stemi.v3.monocytes, features = c(il4.13.v2.Baseline.vs.t24h.mono.targets, il4.13.v3.Baseline.vs.t24h.mono.targets), assays = 'SCT', group.by = 'seurat_clusters')

