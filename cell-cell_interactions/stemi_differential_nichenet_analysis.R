#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_differential_nichenet_analysis.R
# Function: perform cell-cell interaction analysis using Nichenet
############################################################################################################################

####################
# libraries        #
####################

# required for nichenet
library(Seurat)
library(nichenetr)
library(tidyverse)


####################
# Functions        #
####################


calculate_niche_de <- function (seurat_obj, niches, type, assay_oi = "SCT", expression_pct=0.1, logfc_threshold=0.1, test_use='MAST') 
{
  if (type == "sender") {
    sender_vs_sender_tbl = NULL
    for (niche_n in seq(length(niches))) {
      niche = niches[[niche_n]]
      senders_niche = niche$sender %>% unlist() %>% unique()
      senders_other = niches %>% sapply(function(niche) {
        niche$sender
      }) %>% unlist() %>% unique() %>% setdiff(senders_niche)
      for (sender_niche_n in seq(length(senders_niche))) {
        sender_oi = senders_niche[sender_niche_n]
        for (sender_other_n in seq(length(senders_other))) {
          sender_other_niche_oi = senders_other[sender_other_n]
          if (is.null(sender_vs_sender_tbl)) {
            sender_vs_sender_tbl = tibble(sender = sender_oi, 
                                          sender_other_niche = sender_other_niche_oi)
          }
          else {
            if (nrow(sender_vs_sender_tbl %>% filter(sender == 
                                                     sender_other_niche_oi & sender_other_niche == 
                                                     sender_oi)) != 1) {
              sender_vs_sender_tbl = bind_rows(sender_vs_sender_tbl, 
                                               tibble(sender = sender_oi, sender_other_niche = sender_other_niche_oi))
            }
          }
        }
      }
    }
    DE_sender = niches %>% lapply(function(niche, seurat_obj) {
      senders_niche = niche$sender %>% unlist() %>% unique()
      senders_other = niches %>% sapply(function(niche) {
        niche$sender
      }) %>% unlist() %>% unique() %>% setdiff(senders_niche)
      DE_sender = senders_niche %>% lapply(function(sender_oi, 
                                                    seurat_obj, senders_other) {
        print(paste0("Calculate Sender DE between: ", 
                     sender_oi, " and ", senders_other))
        DE_subtable = senders_other %>% lapply(function(sender_other_niche_oi, 
                                                        seurat_obj, sender_oi) {
          if (nrow(sender_vs_sender_tbl %>% filter(sender == 
                                                   sender_other_niche_oi & sender_other_niche == 
                                                   sender_oi)) == 1) {
            DE_sender_oi = NULL
          }
          else {
            DE_sender_oi = FindMarkers(object = seurat_obj, 
                                       ident.1 = sender_oi, ident.2 = sender_other_niche_oi, 
                                       min.pct = expression_pct, logfc.threshold = logfc_threshold, only.pos = FALSE, 
                                       assay = assay_oi, test.use = test_use) %>% rownames_to_column("gene") %>% 
              as_tibble()
            SeuratV4 = c("avg_log2FC") %in% colnames(DE_sender_oi)
            if (SeuratV4 == FALSE) {
              DE_sender_oi = DE_sender_oi %>% dplyr::rename(avg_log2FC = avg_logFC)
            }
            DE_sender_oi = DE_sender_oi %>% mutate(sender = sender_oi, 
                                                   sender_other_niche = sender_other_niche_oi) %>% 
              arrange(-avg_log2FC)
          }
        }, seurat_obj, sender_oi) %>% bind_rows()
      }, seurat_obj, senders_other) %>% bind_rows()
    }, seurat_obj) %>% bind_rows()
    SeuratV4 = c("avg_log2FC") %in% colnames(DE_sender)
    if (SeuratV4 == FALSE) {
      DE_sender = DE_sender %>% dplyr::rename(avg_log2FC = avg_logFC)
    }
    DE_sender_reverse = DE_sender %>% mutate(avg_log2FC = avg_log2FC * 
                                               -1) %>% rename(pct.1_old = pct.1, pct.2_old = pct.2, 
                                                              sender_old = sender, sender_other_niche_old = sender_other_niche) %>% 
      rename(pct.1 = pct.2_old, pct.2 = pct.1_old, sender = sender_other_niche_old, 
             sender_other_niche = sender_old) %>% select(gene, 
                                                         p_val, avg_log2FC, pct.1, pct.2, p_val_adj, sender, 
                                                         sender_other_niche)
    DE_sender = bind_rows(DE_sender, DE_sender_reverse) %>% 
      distinct()
    return(DE_sender)
  }
  if (type == "receiver") {
    receiver_vs_receiver_tbl = NULL
    for (niche_n in seq(length(niches))) {
      niche = niches[[niche_n]]
      receivers_niche = niche$receiver %>% unlist() %>% 
        unique()
      receivers_other = niches %>% sapply(function(niche) {
        niche$receiver
      }) %>% unlist() %>% unique() %>% setdiff(receivers_niche)
      for (receiver_niche_n in seq(length(receivers_niche))) {
        receiver_oi = receivers_niche[receiver_niche_n]
        for (receiver_other_n in seq(length(receivers_other))) {
          receiver_other_niche_oi = receivers_other[receiver_other_n]
          if (is.null(receiver_vs_receiver_tbl)) {
            receiver_vs_receiver_tbl = tibble(receiver = receiver_oi, 
                                              receiver_other_niche = receiver_other_niche_oi)
          }
          else {
            if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == 
                                                         receiver_other_niche_oi & receiver_other_niche == 
                                                         receiver_oi)) != 1) {
              receiver_vs_receiver_tbl = bind_rows(receiver_vs_receiver_tbl, 
                                                   tibble(receiver = receiver_oi, receiver_other_niche = receiver_other_niche_oi))
            }
          }
        }
      }
    }
    print(receiver_vs_receiver_tbl)
    DE_receiver = niches %>% lapply(function(niche, seurat_obj) {
      receivers_niche = niche$receiver %>% unlist() %>% 
        unique()
      receivers_other = niches %>% sapply(function(niche) {
        niche$receiver
      }) %>% unlist() %>% unique() %>% setdiff(receivers_niche)
      DE_receiver = receivers_niche %>% lapply(function(receiver_oi, 
                                                        seurat_obj, receivers_other) {
        print(paste0("Calculate receiver DE between: ", 
                     receiver_oi, " and ", receivers_other))
        DE_subtable = receivers_other %>% lapply(function(receiver_other_niche_oi, 
                                                          seurat_obj, receiver_oi) {
          if (nrow(receiver_vs_receiver_tbl %>% filter(receiver == 
                                                       receiver_other_niche_oi & receiver_other_niche == 
                                                       receiver_oi)) == 1) {
            DE_receiver_oi = NULL
          }
          else {
            DE_receiver_oi = FindMarkers(object = seurat_obj, 
                                         ident.1 = receiver_oi, ident.2 = receiver_other_niche_oi, 
                                         min.pct = expression_pct, logfc.threshold = logfc_threshold, only.pos = FALSE, 
                                         assay = assay_oi, test.use = test_use) %>% rownames_to_column("gene") %>% 
              as_tibble()
            SeuratV4 = c("avg_log2FC") %in% colnames(DE_receiver_oi)
            if (SeuratV4 == FALSE) {
              DE_receiver_oi = DE_receiver_oi %>% dplyr::rename(avg_log2FC = avg_logFC)
            }
            DE_receiver_oi = DE_receiver_oi %>% mutate(receiver = receiver_oi, 
                                                       receiver_other_niche = receiver_other_niche_oi) %>% 
              arrange(-avg_log2FC)
          }
        }, seurat_obj, receiver_oi) %>% bind_rows()
      }, seurat_obj, receivers_other) %>% bind_rows()
    }, seurat_obj) %>% bind_rows()
    SeuratV4 = c("avg_log2FC") %in% colnames(DE_receiver)
    if (SeuratV4 == FALSE) {
      DE_receiver = DE_receiver %>% dplyr::rename(avg_log2FC = avg_logFC)
    }
    DE_receiver_reverse = DE_receiver %>% mutate(avg_log2FC = avg_log2FC * 
                                                   -1) %>% rename(pct.1_old = pct.1, pct.2_old = pct.2, 
                                                                  receiver_old = receiver, receiver_other_niche_old = receiver_other_niche) %>% 
      rename(pct.1 = pct.2_old, pct.2 = pct.1_old, receiver = receiver_other_niche_old, 
             receiver_other_niche = receiver_old) %>% select(gene, 
                                                             p_val, avg_log2FC, pct.1, pct.2, p_val_adj, receiver, 
                                                             receiver_other_niche)
    DE_receiver = bind_rows(DE_receiver, DE_receiver_reverse) %>% 
      distinct()
    return(DE_receiver)
  }
}



do_differential_nichenet <- function(seurat_obj, celltype_column, condition_column, condition_1, condition_2, cell_types_of_interest=NULL, assay_oi='SCT', expression_pct=0.1, logfc_threshold=0.1, test_use='MAST', specificity_score_LR_pairs="min_lfc", lfc_cutoff=0.15, top_n_target=250) {
  # subset the seurat object
  seurat_obj <- seurat_obj[, seurat_obj@meta.data[[condition_column]] %in% c(condition_1, condition_2)]
  # add a new column to have both the cell type and the condition
  seurat_obj@meta.data[['ct_cond']] <- paste(seurat_obj@meta.data[[celltype_column]], seurat_obj@meta.data[[condition_column]], sep = '_')
  # set that as the idents
  Idents(seurat_obj) <- 'ct_cond'
  # get the cell types of interest
  if (is.null(cell_types_of_interest)) {
    # setting to all if not supplied
    cell_types_of_interest <- unique(seurat_obj@meta.data[[celltype_column]])
  }
  # get the niches list
  niches <- list()
  niches[[condition_1]] <- list(
    'sender' = paste(cell_types_of_interest, condition_1, sep = '_'),
    'receiver' = paste(cell_types_of_interest, condition_1, sep = '_')
  )
  niches[[condition_2]] <- list(
    'sender' = paste(cell_types_of_interest, condition_2, sep = '_'),
    'receiver' = paste(cell_types_of_interest, condition_2, sep = '_')
  )
  # convert symbols
  seurat_obj = alias_to_symbol_seurat(seurat_obj, organism = "human")
  # get the sending niches
  DE_sender = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$ligand %>% unique()), niches = niches, type = "sender", assay_oi = assay_oi, expression_pct=expression_pct, logfc_threshold=logfc_threshold, test_use=test_use)
  # get the receiving niches
  DE_receiver = calculate_niche_de(seurat_obj = seurat_obj %>% subset(features = lr_network$receptor %>% unique()), niches = niches, type = "receiver", assay_oi = assay_oi, expression_pct=expression_pct, logfc_threshold=logfc_threshold, test_use=test_use)
  # remove infinite LFC, where the expression was zero in either conditions (relevant but computationally difficult to work with)
  DE_sender = DE_sender %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  DE_receiver = DE_receiver %>% mutate(avg_log2FC = ifelse(avg_log2FC == Inf, max(avg_log2FC[is.finite(avg_log2FC)]), ifelse(avg_log2FC == -Inf, min(avg_log2FC[is.finite(avg_log2FC)]), avg_log2FC)))
  # process DE output
  DE_sender_processed = process_niche_de(DE_table = DE_sender, niches = niches, expression_pct = expression_pct, type = "sender")
  DE_receiver_processed = process_niche_de(DE_table = DE_receiver, niches = niches, expression_pct = expression_pct, type = "receiver")
  # combine the sender and receiver
  DE_sender_receiver = combine_sender_receiver_de(DE_sender_processed, DE_receiver_processed, lr_network, specificity_score = specificity_score_LR_pairs)
  # mock the spatial info, we don't have this, and spatial info doesn't say much in blood
  spatial_info = tibble(celltype_region_oi = NA, celltype_other_region = NA) %>% mutate(niche =  niches %>% names() %>% head(1), celltype_type = "sender")
  # calculate  niches
  DE_receiver_targets = calculate_niche_de_targets(seurat_obj = seurat_obj, niches = niches, lfc_cutoff = lfc_cutoff, expression_pct = expression_pct, assay_oi = assay_oi) 
  DE_receiver_processed_targets = process_receiver_target_de(DE_receiver_targets = DE_receiver_targets, niches = niches, expression_pct = expression_pct, specificity_score = 'min_lfc')
  
  print('fetching gene sets')
  # get background set
  background = DE_receiver_processed_targets  %>% pull(target) %>% unique()
  # get the gene sets
  geneset_niche1 = DE_receiver_processed_targets %>% filter(receiver == niches[[1]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  geneset_niche2 = DE_receiver_processed_targets %>% filter(receiver == niches[[2]]$receiver & target_score >= lfc_cutoff & target_significant == 1 & target_present == 1) %>% pull(target) %>% unique()
  # filtering for unmappable names
  geneset_niche1 %>% setdiff(rownames(ligand_target_matrix))
  geneset_niche2 %>% setdiff(rownames(ligand_target_matrix))
  # create the geneset lists, 
  niche_geneset_list = list(
    condition_1 = list(
      "receiver" = niches[[1]]$receiver,
      "geneset" = geneset_niche1,
      "background" = background),
    condition_2 = list(
      "receiver" = niches[[2]]$receiver,
      "geneset" = geneset_niche2 ,
      "background" = background)
  )
  print(niches[[1]]$receiver)
  print('fetching ligand activity targets')
  # now get the ligand-target activities
  ligand_activities_targets = get_ligand_activities_targets(niche_geneset_list = niche_geneset_list, ligand_target_matrix = ligand_target_matrix, top_n_target = top_n_target)
  
  print('fetching features of interest')
  # get the features we aare going to use
  features_oi = union(lr_network$ligand, lr_network$receptor) %>% union(ligand_activities_targets$target) %>% setdiff(NA)
  # use the Seurat dotplot function to get the expression across the idents
  dotplot = suppressWarnings(Seurat::DotPlot(seurat_obj %>% subset(idents = niches %>% unlist() %>% unique()), features = features_oi, assay = assay_oi))
  exprs_tbl = dotplot$data %>% as_tibble()
  exprs_tbl = exprs_tbl %>% rename(celltype = id, gene = features.plot, expression = avg.exp, expression_scaled = avg.exp.scaled, fraction = pct.exp) %>%
    mutate(fraction = fraction/100) %>% as_tibble() %>% select(celltype, gene, expression, expression_scaled, fraction) %>% distinct() %>% arrange(gene) %>% mutate(gene = as.character(gene))
  exprs_tbl_ligand = exprs_tbl %>% filter(gene %in% lr_network$ligand) %>% rename(sender = celltype, ligand = gene, ligand_expression = expression, ligand_expression_scaled = expression_scaled, ligand_fraction = fraction) 
  exprs_tbl_receptor = exprs_tbl %>% filter(gene %in% lr_network$receptor) %>% rename(receiver = celltype, receptor = gene, receptor_expression = expression, receptor_expression_scaled = expression_scaled, receptor_fraction = fraction)
  exprs_tbl_target = exprs_tbl %>% filter(gene %in% ligand_activities_targets$target) %>% rename(receiver = celltype, target = gene, target_expression = expression, target_expression_scaled = expression_scaled, target_fraction = fraction)
  exprs_tbl_ligand = exprs_tbl_ligand %>%  mutate(scaled_ligand_expression_scaled = scale_quantile_adapted(ligand_expression_scaled)) %>% mutate(ligand_fraction_adapted = ligand_fraction) %>% mutate_cond(ligand_fraction >= expression_pct, ligand_fraction_adapted = expression_pct)  %>% mutate(scaled_ligand_fraction_adapted = scale_quantile_adapted(ligand_fraction_adapted))
  exprs_tbl_receptor = exprs_tbl_receptor %>% mutate(scaled_receptor_expression_scaled = scale_quantile_adapted(receptor_expression_scaled))  %>% mutate(receptor_fraction_adapted = receptor_fraction) %>% mutate_cond(receptor_fraction >= expression_pct, receptor_fraction_adapted = expression_pct)  %>% mutate(scaled_receptor_fraction_adapted = scale_quantile_adapted(receptor_fraction_adapted))
  
  print('scoreing L-R interactions')
  # score ligand-receptor interactions based on expression strength of the receptor
  exprs_sender_receiver = lr_network %>% 
    inner_join(exprs_tbl_ligand, by = c("ligand")) %>% 
    inner_join(exprs_tbl_receptor, by = c("receptor")) %>% inner_join(DE_sender_receiver %>% distinct(niche, sender, receiver))
  ligand_scaled_receptor_expression_fraction_df = exprs_sender_receiver %>% group_by(ligand, receiver) %>% mutate(rank_receptor_expression = dense_rank(receptor_expression), rank_receptor_fraction  = dense_rank(receptor_fraction)) %>% mutate(ligand_scaled_receptor_expression_fraction = 0.5*( (rank_receptor_fraction / max(rank_receptor_fraction)) + ((rank_receptor_expression / max(rank_receptor_expression))) ) )  %>% distinct(ligand, receptor, receiver, ligand_scaled_receptor_expression_fraction) %>% distinct() %>% ungroup()
  
  # weight different variables for the prioritization
  prioritizing_weights = c("scaled_ligand_score" = 5,
                           "scaled_ligand_expression_scaled" = 1,
                           "ligand_fraction" = 1,
                           "scaled_ligand_score_spatial" = 2, 
                           "scaled_receptor_score" = 0.5,
                           "scaled_receptor_expression_scaled" = 0.5,
                           "receptor_fraction" = 1, 
                           "ligand_scaled_receptor_expression_fraction" = 1,
                           "scaled_receptor_score_spatial" = 0,
                           "scaled_activity" = 0,
                           "scaled_activity_normalized" = 1)
  
  # get output and priorization
  output = list(DE_sender_receiver = DE_sender_receiver, ligand_scaled_receptor_expression_fraction_df = ligand_scaled_receptor_expression_fraction_df, sender_spatial_DE_processed = sender_spatial_DE_processed, receiver_spatial_DE_processed = receiver_spatial_DE_processed,
                ligand_activities_targets = ligand_activities_targets, DE_receiver_processed_targets = DE_receiver_processed_targets, exprs_tbl_ligand = exprs_tbl_ligand,  exprs_tbl_receptor = exprs_tbl_receptor, exprs_tbl_target = exprs_tbl_target)
  prioritization_tables = get_prioritization_tables(output, prioritizing_weights)
  
  return(prioritization_tables)
}


####################
# Main Code        #
####################

# location of the nichenet database files
nichenet_database_locs_base <- '/groups/umcg-franke-scrna/tmp01/external_datasets/nichenet/'
ligand_target_matrix_loc <- paste(nichenet_database_locs_base, 'ligand_target_matrix_nsga2r_final.rds', sep = '')
lr_network_loc <- paste(nichenet_database_locs_base, 'lr_network_human_21122021.rds', sep = '')

# location of the object
objects_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/'
cardio_integrated_loc <- paste(objects_loc, 'cardio.integrated.20210301.rds', sep = '')

# read the object
cardio_integrated <- readRDS(cardio_integrated_loc)
cardio_integrated <- UpdateSeuratObject(cardio_integrated)
DefaultAssay(cardio_integrated) <- 'SCT'

# read the databases
ligand_target_matrix <- readRDS(ligand_target_matrix_loc)
lr_network_loc <- readRDS(lr_network_loc)
lr_network = lr_network %>% dplyr::rename(ligand = from, receptor = to) %>% distinct(ligand, receptor)

# do prioritization of t0 vs t8w
prio_t0_t8w <- do_differential_nichenet(cardio_integrated, 'cell_type_lowerres', 'timepoint.final', 't8w', 'Baseline', test_use='wilcox', cell_types_of_interest = c('B', 'monocyte', 'NK'))
