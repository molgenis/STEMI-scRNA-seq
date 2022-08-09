#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_plot_kegg_ctc_pathways.R
# Function:
############################################################################################################################

####################
# libraries        #
####################


library(organism, character.only = TRUE)
library(enrichplot)
library(pathview)
library(ggnewscale)
library(clusterProfiler)

####################
# Functions        #
####################


extract_significant_ligands_and_downstream_genes <- function(nichenet_output, significance_cutoff=0.05, matrix='ligand_target'){
  # make a list per receiver
  interactions <- list()
  # check each receiver
  for(receiver in names(nichenet_output)){
    # make a list per sender
    interactions_senders <- list()
    # check against each sender
    for(sender in names(nichenet_output[[receiver]])){
      # init lists of ligands and downstream genes
      ligands_significant <- c()
      downstream_genes_significant <- c()
      # extract the ligand-downstream gene matrix
      if(matrix %in% names(nichenet_output[[receiver]][[sender]])){
        ligand_gene_matrix <- nichenet_output[[receiver]][[sender]][[matrix]]
        if(nrow(ligand_gene_matrix) > 0){
          # extract the ligands
          ligands <- colnames(ligand_gene_matrix)
          # extract the downstream genes
          downstream_genes <- rownames(ligand_gene_matrix)
          # check each combination
          for(ligand in ligands){
            for(downstream_gene in downstream_genes){
               # get the value
              interaction <- ligand_gene_matrix[downstream_gene, ligand]
              # if larger than the cutoff value (0.05 is nichenet default)
              if(interaction > significance_cutoff){
                # add both the ligand and the downstream gene
                ligands_significant <- c(ligands_significant, ligand)
                downstream_genes_significant <- c(downstream_genes_significant, downstream_gene)
              }
            }
          }
        }
      }
      # put genes is list
      ligand_targets <- list('ligands' = ligands_significant, 'targets' = downstream_genes_significant)
      # save list for the sender
      interactions_senders[[sender]] <- ligand_targets
    }
    # save list for receiver
    interactions[[receiver]] <- interactions_senders
  }
  return(interactions)
}

get_unique_ligands_targets <- function(significant_ligands_targets, ligand_column='ligands', target_column='targets'){
  # init lists of ligands and downstream genes
  ligands_significant <- c()
  downstream_genes_significant <- c()
  # check each receiver
  for(receiver in names(significant_ligands_targets)){
    # check against each sender
    for(sender in names(significant_ligands_targets[[receiver]])){
      # extract these ligands
      ligands <- significant_ligands_targets[[receiver]][[sender]][[ligand_column]]
      # extract the downstream genes
      downstream_genes <- significant_ligands_targets[[receiver]][[sender]][[target_column]]
      # add to the total list
      ligands_significant <- c(ligands_significant, ligands)
      downstream_genes_significant <- c(downstream_genes_significant, downstream_genes)
    }
  }
  return(list('ligands' = unique(ligands_significant), 'targets' = unique(downstream_genes_significant)))
}

get_lfs_for_genes <- function(mast_output_loc, condition_combination, genes, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), lfc_column='metafc'){
  lfc_dataframe <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # paste the output path together
    full_output_path <- paste(mast_output_loc, cell_type, condition_combination, '.tsv', sep = '')
    # see if the file exists
    if(file.exists(full_output_path)){
      # read the file
      de_output <- read.table(full_output_path, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
      # get the lfs
      lfcs <- de_output[match(genes, rownames(de_output)), lfc_column]
      # turn into a dataframe
      lfcs_column <- data.frame(lfcs=lfcs)
      colnames(lfcs_column) <- cell_type
      # add to the dataframe
      if(is.null(lfc_dataframe)){
        rownames(lfcs_column) <- genes
        lfc_dataframe <- lfcs_column
      }
      else{
        lfc_dataframe <- cbind(lfc_dataframe, lfcs_column)
      }
    }
  }
  return(lfc_dataframe)
}


get_absolute_max <- function(vector){
  abs_max <- vector[which.max(abs(vector))]
  return(abs_max)
}

####################
# Main Code        #
####################


# set the organism
kegg.organism <- "hsa"
organism <- 'org.Hs.eg.db'
# mast output location
mast_output_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_20200707/rna/'
# read Baseline vs 24h
v2_Baseline_vs_t24h_perct_omni_unweighted_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v2_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted.rds'
v2_Baseline_vs_t24h_perct_omni_unweighted <- readRDS(v2_Baseline_vs_t24h_perct_omni_unweighted_loc)
v2_Baseline_vs_t24h_perct_omni_unweighted_genes <- extract_significant_ligands_and_downstream_genes(v2_Baseline_vs_t24h_perct_omni_unweighted)
v2_Baseline_vs_t24h_perct_omni_unweighted_genes_all <- get_unique_ligands_targets(v2_Baseline_vs_t24h_perct_omni_unweighted_genes)
v3_Baseline_vs_t24h_perct_omni_unweighted_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v3_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted.rds'
v3_Baseline_vs_t24h_perct_omni_unweighted <- readRDS(v3_Baseline_vs_t24h_perct_omni_unweighted_loc)
v3_Baseline_vs_t24h_perct_omni_unweighted_genes <- extract_significant_ligands_and_downstream_genes(v3_Baseline_vs_t24h_perct_omni_unweighted)
v3_Baseline_vs_t24h_perct_omni_unweighted_genes_all <- get_unique_ligands_targets(v3_Baseline_vs_t24h_perct_omni_unweighted_genes)
# now also the receptors
v2_Baseline_vs_t24h_perct_omni_unweighted_receptors <- extract_significant_ligands_and_downstream_genes(v2_Baseline_vs_t24h_perct_omni_unweighted, matrix = 'ligand_receptor_network')
v2_Baseline_vs_t24h_perct_omni_unweighted_receptors_all <- get_unique_ligands_targets(v2_Baseline_vs_t24h_perct_omni_unweighted_receptors)
v3_Baseline_vs_t24h_perct_omni_unweighted_receptors <- extract_significant_ligands_and_downstream_genes(v3_Baseline_vs_t24h_perct_omni_unweighted, matrix = 'ligand_receptor_network')
v3_Baseline_vs_t24h_perct_omni_unweighted_receptors_all <- get_unique_ligands_targets(v3_Baseline_vs_t24h_perct_omni_unweighted_receptors)

# merge the version chemistry
Baseline_vs_t24h_perct_omni_unweighted_ligands <- unique(c(v2_Baseline_vs_t24h_perct_omni_unweighted_genes_all[['ligands']],
                                                           v3_Baseline_vs_t24h_perct_omni_unweighted_genes_all[['ligands']],
                                                           v2_Baseline_vs_t24h_perct_omni_unweighted_receptors_all[['ligands']],
                                                           v3_Baseline_vs_t24h_perct_omni_unweighted_receptors_all[['ligands']]))
Baseline_vs_t24h_perct_omni_unweighted_targets <- unique(c(v2_Baseline_vs_t24h_perct_omni_unweighted_genes_all[['targets']],
                                                           v3_Baseline_vs_t24h_perct_omni_unweighted_genes_all[['targets']],
                                                           v2_Baseline_vs_t24h_perct_omni_unweighted_receptors_all[['targets']],
                                                           v3_Baseline_vs_t24h_perct_omni_unweighted_receptors_all[['targets']]))
# combine the ligands and downstream genes
Baseline_vs_t24h_perct_omni_unweighted_ligands_targets <- c(Baseline_vs_t24h_perct_omni_unweighted_ligands, Baseline_vs_t24h_perct_omni_unweighted_targets)
# get the entrez IDs for the gene symbols
Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping <- bitr(Baseline_vs_t24h_perct_omni_unweighted_ligands_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# get the IDs which are not duplicated, and are not NA
Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping <- Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping[!is.na(Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping[['ENTREZID']]), ]
# now get the LFCs
Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_lfcs <- get_lfs_for_genes(mast_output_loc, 'Baselinet24h', Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping[['SYMBOL']])
# set to zero if empty
Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_lfcs[is.na(Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_lfcs)] <- 0
# get the max absolute LFC
Baseline_vs_t24h_perct_omni_unweighted_max_lfcs <- apply(Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_lfcs, 1, get_absolute_max)
# extract the values
Baseline_vs_t24h_perct_omni_unweighted_genes_lfcs <- names(Baseline_vs_t24h_perct_omni_unweighted_max_lfcs)
Baseline_vs_t24h_perct_omni_unweighted_values_lfcs <- as.vector(unlist(Baseline_vs_t24h_perct_omni_unweighted_max_lfcs))
Baseline_vs_t24h_perct_omni_unweighted_genes_lfcs <- Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping[match(Baseline_vs_t24h_perct_omni_unweighted_genes_lfcs, Baseline_vs_t24h_perct_omni_unweighted_ligands_targets_mapping[['SYMBOL']]), 'ENTREZID']
names(Baseline_vs_t24h_perct_omni_unweighted_values_lfcs) <- Baseline_vs_t24h_perct_omni_unweighted_genes_lfcs
# view the pathways
wd_prev <- getwd()
dir.create(paste(wd_prev, '/Baselinet24h/', sep = ''))
setwd(paste(wd_prev, '/Baselinet24h/', sep = ''))
pathview(gene.data=Baseline_vs_t24h_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04060', species = kegg.organism)
pathview(gene.data=Baseline_vs_t24h_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04061', species = kegg.organism)
pathview(gene.data=Baseline_vs_t24h_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04657', species = kegg.organism)
pathview(gene.data=Baseline_vs_t24h_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04659', species = kegg.organism)
setwd(wd_prev)

# read Baseline vs 24h
v2_Baseline_vs_t8w_perct_omni_unweighted_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v2_Baseline_vs_t8w_nichenet_onlymajor_perct_omni_unweighted.rds'
v2_Baseline_vs_t8w_perct_omni_unweighted <- readRDS(v2_Baseline_vs_t8w_perct_omni_unweighted_loc)
v2_Baseline_vs_t8w_perct_omni_unweighted_genes <- extract_significant_ligands_and_downstream_genes(v2_Baseline_vs_t8w_perct_omni_unweighted)
v2_Baseline_vs_t8w_perct_omni_unweighted_genes_all <- get_unique_ligands_targets(v2_Baseline_vs_t8w_perct_omni_unweighted_genes)
v3_Baseline_vs_t8w_perct_omni_unweighted_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/objects/v3_Baseline_vs_t8w_nichenet_onlymajor_perct_omni_unweighted.rds'
v3_Baseline_vs_t8w_perct_omni_unweighted <- readRDS(v3_Baseline_vs_t8w_perct_omni_unweighted_loc)
v3_Baseline_vs_t8w_perct_omni_unweighted_genes <- extract_significant_ligands_and_downstream_genes(v3_Baseline_vs_t8w_perct_omni_unweighted)
v3_Baseline_vs_t8w_perct_omni_unweighted_genes_all <- get_unique_ligands_targets(v3_Baseline_vs_t8w_perct_omni_unweighted_genes)
# now also the receptors
v2_Baseline_vs_t8w_perct_omni_unweighted_receptors <- extract_significant_ligands_and_downstream_genes(v2_Baseline_vs_t8w_perct_omni_unweighted, matrix = 'ligand_receptor_network')
v2_Baseline_vs_t8w_perct_omni_unweighted_receptors_all <- get_unique_ligands_targets(v2_Baseline_vs_t8w_perct_omni_unweighted_receptors)
v3_Baseline_vs_t8w_perct_omni_unweighted_receptors <- extract_significant_ligands_and_downstream_genes(v3_Baseline_vs_t8w_perct_omni_unweighted, matrix = 'ligand_receptor_network')
v3_Baseline_vs_t8w_perct_omni_unweighted_receptors_all <- get_unique_ligands_targets(v3_Baseline_vs_t8w_perct_omni_unweighted_receptors)

# merge the version chemistry
Baseline_vs_t8w_perct_omni_unweighted_ligands <- unique(c(v2_Baseline_vs_t8w_perct_omni_unweighted_genes_all[['ligands']],
                                                           v3_Baseline_vs_t8w_perct_omni_unweighted_genes_all[['ligands']],
                                                           v2_Baseline_vs_t8w_perct_omni_unweighted_receptors_all[['ligands']],
                                                           v3_Baseline_vs_t8w_perct_omni_unweighted_receptors_all[['ligands']]))
Baseline_vs_t8w_perct_omni_unweighted_targets <- unique(c(v2_Baseline_vs_t8w_perct_omni_unweighted_genes_all[['targets']],
                                                           v3_Baseline_vs_t8w_perct_omni_unweighted_genes_all[['targets']],
                                                           v2_Baseline_vs_t8w_perct_omni_unweighted_receptors_all[['targets']],
                                                           v3_Baseline_vs_t8w_perct_omni_unweighted_receptors_all[['targets']]))
# combine the ligands and downstream genes
Baseline_vs_t8w_perct_omni_unweighted_ligands_targets <- c(Baseline_vs_t8w_perct_omni_unweighted_ligands, Baseline_vs_t8w_perct_omni_unweighted_targets)
# get the entrez IDs for the gene symbols
Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping <- bitr(Baseline_vs_t8w_perct_omni_unweighted_ligands_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
# get the IDs which are not duplicated, and are not NA
Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping <- Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping[!is.na(Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping[['ENTREZID']]), ]
# now get the LFCs
Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_lfcs <- get_lfs_for_genes(mast_output_loc, 'Baselinet8w', Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping[['SYMBOL']])
# set to zero if empty
Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_lfcs[is.na(Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_lfcs)] <- 0
# get the max absolute LFC
Baseline_vs_t8w_perct_omni_unweighted_max_lfcs <- apply(Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_lfcs, 1, get_absolute_max)
# extract the values
Baseline_vs_t8w_perct_omni_unweighted_genes_lfcs <- names(Baseline_vs_t8w_perct_omni_unweighted_max_lfcs)
Baseline_vs_t8w_perct_omni_unweighted_values_lfcs <- as.vector(unlist(Baseline_vs_t8w_perct_omni_unweighted_max_lfcs))
Baseline_vs_t8w_perct_omni_unweighted_genes_lfcs <- Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping[match(Baseline_vs_t8w_perct_omni_unweighted_genes_lfcs, Baseline_vs_t8w_perct_omni_unweighted_ligands_targets_mapping[['SYMBOL']]), 'ENTREZID']
names(Baseline_vs_t8w_perct_omni_unweighted_values_lfcs) <- Baseline_vs_t8w_perct_omni_unweighted_genes_lfcs
# view the pathways
wd_prev <- getwd()
dir.create(paste(wd_prev, '/Baselinet8w/', sep = ''))
setwd(paste(wd_prev, '/Baselinet8w/', sep = ''))
pathview(gene.data=Baseline_vs_t8w_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04060', species = kegg.organism)
pathview(gene.data=Baseline_vs_t8w_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04061', species = kegg.organism)
pathview(gene.data=Baseline_vs_t8w_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04657', species = kegg.organism)
pathview(gene.data=Baseline_vs_t8w_perct_omni_unweighted_values_lfcs, pathway.id = 'hsa04659', species = kegg.organism)
setwd(wd_prev)
