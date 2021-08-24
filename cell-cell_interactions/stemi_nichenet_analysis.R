#
# some description of the script
#

#
# libraries
#

library(nichenetr)
library(Seurat)
library(tidyverse)

#
# functions
#


do_nichenet_analysis <- function(seurat_object, receiver, sender_celltypes, condition.column, condition.1, condition.2, lr_network, weighted_networks, ligand_target_matrix, cell_type_column='cell_type_lowerres', pct=0.1, min_avg_log2FC=0.25, top_n=20, max_active_ligand_target_links=200, active_ligand_target_links_cutoff=0.33, only_documented_lr=F, test.use='MAST'){
  # we'll do a couple of analyses, so we'll store the results in a list
  analysis_list <- list()
  # weigh network
  weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
  # if requested, don't use predicted, but only documented LR interactions
  if(only_documented_lr){
    lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
  }
  # init receiver data
  expressed_genes_receiver = get_expressed_genes(receiver, seurat_object, pct = pct)
  background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  # init sender data
  list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_object, pct) # lapply to get the expressed genes of every sender cell type separately here
  expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
  # subset to receiver
  seurat_obj_receiver= seurat_object[, seurat_object@meta.data[[cell_type_column]] == receiver]
  # set condition column as ident for findmarkers
  seurat_obj_receiver = SetIdent(seurat_obj_receiver, value = seurat_obj_receiver[[condition.column]])
  condition_oi = condition.1
  condition_reference = condition.2
  try({
    DE_table_receiver = FindMarkers(object = seurat_obj_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = pct, test.use = test.use) %>% tibble::rownames_to_column("gene")
    # subset to significant DE genes
    geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= min_avg_log2FC) %>% pull(gene)
    geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
    # get the ligands and receptors
    ligands = lr_network %>% pull(from) %>% unique()
    receptors = lr_network %>% pull(to) %>% unique()
    # subset to expressed ones
    expressed_ligands = intersect(ligands,expressed_genes_sender)
    expressed_receptors = intersect(receptors,expressed_genes_receiver)
    # subset to potential ligands by including only the expressed ones
    potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
    # do the ligand activity analysis
    ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
    # sort by correlation
    ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
    # subset to top ligands
    best_upstream_ligands = ligand_activities %>% top_n(top_n, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
    analysis_list[['best_upstream_ligands']] <- best_upstream_ligands
    # infer receptors and top-predicted target genes of top-ranked ligands
    active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = max_active_ligand_target_links) %>% bind_rows() %>% tidyr::drop_na()
    active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = active_ligand_target_links_cutoff)
    # get ligands for which there are links
    order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
    # get the targets for which there are links
    order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
    # create vis-compatible matrix
    rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
    colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
    # add to list
    analysis_list[['ligand_target']] <- vis_ligand_target
    # filter the L-R network by taking the best ligands that have expressed receptors
    lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
    # get the best unique upstream receptors
    best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()
    # filter to the best ligands and best receptors
    lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
    # set the combinations of lr to be one ligand and the receptors as columns
    lr_network_top_df = lr_network_top_df_large %>% pivot_wider(names_from = "from", values_from = "weight",values_fill = 0)
    # ?not completely sure here
    lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
    # cluster weights
    dist_receptors = dist(lr_network_top_matrix, method = "binary")
    hclust_receptors = hclust(dist_receptors, method = "ward.D2")
    # extract the receptor order
    order_receptors = hclust_receptors$labels[hclust_receptors$order]
    # cluster again
    dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary") # same as before actually
    hclust_ligands = hclust(dist_ligands, method = "ward.D2")
    # extract order again
    order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]
    # subset to LR present in top matrix
    order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
    order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))
    # create vis-compatible matrix again
    vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
    rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
    colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()
    # add to list
    analysis_list[['ligand_receptor_network']] <- vis_ligand_receptor_network
    # add full DE data
    DE_table_all <- c(sender_celltypes)  %>% lapply(get_lfc_celltype, seurat_obj = seurat_object, condition_colname = condition.column, condition_oi = condition.2, condition_reference = condition.1, expression_pct = pct,test.use = test.use) %>% reduce(full_join)
    DE_table_all[is.na(DE_table_all)] = 0
    # Combine ligand activities with DE information
    ligand_activities_de = ligand_activities %>% select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
    ligand_activities_de[is.na(ligand_activities_de)] = 0
    # make LFC heatmap
    lfc_matrix = ligand_activities_de  %>% select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
    rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()
    # order by the LFC
    order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
    vis_ligand_lfc = lfc_matrix[order_ligands,]
    colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()
    # add to list
    analysis_list[['ligand_lfc']] <- vis_ligand_lfc
  })
  return(analysis_list)
}


do_nichenet_analysis_per_celltype <- function(seurat_object, condition.column, condition.1, condition.2, lr_network, weighted_networks, ligand_target_matrix, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column='cell_type_lowerres', pct=0.1, min_avg_log2FC=0.25, top_n=20, max_active_ligand_target_links=200, active_ligand_target_links_cutoff=0.33, only_documented_lr=F){
  # confine to cell types that we have in our dataset
  cell_types_to_use <- intersect(cell_types, unique(seurat_object@meta.data[[cell_type_column]]))
  # set cell type as the ident
  Idents(seurat_object) <- cell_type_column
  # and add as separate column
  seurat_object@meta.data[['celltype']] <- seurat_object@meta.data[[cell_type_column]]
  # save results per cell type
  analysis_per_cell_type <- list()
  # check each cell type against the others
  for(cell_type in cell_types_to_use){
    # get the other cell types
    cell_types_other <- setdiff(cell_types_to_use, cell_type)
    # do the analysis
    results <- do_nichenet_analysis(seurat_object, cell_type, cell_types_other, condition.column, condition.1, condition.2, lr_network, weighted_networks, ligand_target_matrix, cell_type_column=cell_type_column, pct=pct, min_avg_log2FC=min_avg_log2FC, top_n=top_n, max_active_ligand_target_links=max_active_ligand_target_links, active_ligand_target_links_cutoff=active_ligand_target_links_cutoff, only_documented_lr=only_documented_lr)
    # add result to list
    analysis_per_cell_type[[cell_type]] <- results
  }
  return(analysis_per_cell_type)
}


do_nichenet_analysis_versus_each_celltype <- function(seurat_object, condition.column, condition.1, condition.2, lr_network, weighted_networks, ligand_target_matrix, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column='cell_type_lowerres', pct=0.1, min_avg_log2FC=0.25, top_n=20, max_active_ligand_target_links=200, active_ligand_target_links_cutoff=0.33, only_documented_lr=F){
  # confine to cell types that we have in our dataset
  cell_types_to_use <- intersect(cell_types, unique(seurat_object@meta.data[[cell_type_column]]))
  # set cell type as the ident
  Idents(seurat_object) <- cell_type_column
  # and add as separate column
  seurat_object@meta.data[['celltype']] <- seurat_object@meta.data[[cell_type_column]]
  # store per cell type
  result_per_cell_type <- list()
  # check each cell type
  for(cell_type1 in cell_types_to_use){
    # store for this cell type
    result_against_cell_type <- list()
    # against each other cell type
    for(cell_type2 in setdiff(cell_types_to_use, cell_type1)){
      # subset to those cell types
      seurat_object_cell_types <- seurat_object[, seurat_object@meta.data$celltype %in% c(cell_type1, cell_type2)]
      # do the analysis for these two
      results <- do_nichenet_analysis(seurat_object=seurat_object_cell_types, cell_type=cell_type1, cell_types_other=c(cell_type2), condition.column=condition.column, condition.1=condition.1, condition.2=condition.2, lr_network=lr_network, weighted_networks=weighted_networks, ligand_target_matrix=ligand_target_matrix, cell_type_column=cell_type_column, pct=pct, min_avg_log2FC=min_avg_log2FC, top_n=top_n, max_active_ligand_target_links=max_active_ligand_target_links, active_ligand_target_links_cutoff=active_ligand_target_links_cutoff, only_documented_lr=only_documented_lr)
      # store the result
      result_against_cell_type[[cell_type2]] <- results
    }
    # store every result
    result_per_cell_type[[cell_type]] <- result_against_cell_type
  }
  return(result_per_cell_type)
}


nichenet_output_to_plot <- function(output_list){
  # store per celltype the resulting plot
  combined_plots <- list()
  # we have multiple cell types
  for(cell_type in names(output_list)){
    # grab that result
    output_list_celltype <- output_list[[cell_type]]

    # create a list to store ggplot2 objects
    plots <- list()

    if((length(output_list_celltype) > 0) & ('ligand_target' %in% names(output_list_celltype))){
      # fetch data from list
      vis_ligand_target <- output_list_celltype[['ligand_target']]
      if(nrow(vis_ligand_target) > 0){
        plots[['ligand_target']] <- vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
      }
    }
    if(length(output_list_celltype) > 0 & 'ligand_receptor_network' %in% names(output_list_celltype)){
      vis_ligand_receptor_network <- output_list_celltype[['ligand_receptor_network']]
      if(nrow(vis_ligand_receptor_network) > 0){
        plots[['ligand_receptor_network']] <- vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
      }
    }
    if(length(output_list_celltype) > 0 & 'ligand_lfc' %in% names(output_list_celltype)){
      vis_ligand_lfc <- output_list_celltype[['ligand_lfc']]
      if(nrow(vis_ligand_lfc) > 0){
        plots[['ligand_lfc']] <- vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
      }
    }
    # combined the plots of the cell type
    if(length(plots) > 0){
      plot_combined <- cowplot::plot_grid(plotlist = plots)
    }
    # add resulting plot to list
    combined_plots[[cell_type]] <- plot_combined
  }
  return(combined_plots)
}

perct_output_to_plot <- function(output_list){
  # store per celltype1 the resulting plot
  combined_plots_ct1 <- list()
  # we have multiple cell types
  for(cell_type1 in names(output_list)){
    # grab that result
    output_list_celltype1 <- output_list[[cell_type1]]
    # store per celltype the resulting plot
    combined_plots_ct2 <- list()
    #get the cell type the interactions are with
    for(cell_type2 in names(output_list_celltype1)){
      output_list_celltype2 <- output_list_celltype1[[cell_type2]]
      # create a list to store ggplot2 objects
      plots <- list()
      
      if((length(output_list_celltype2) > 0) & ('ligand_target' %in% names(output_list_celltype2))){
        # fetch data from list
        vis_ligand_target <- output_list_celltype2[['ligand_target']]
        if(nrow(vis_ligand_target) > 0){
          plots[['ligand_target']] <- vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
        }
      }
      if(length(output_list_celltype2) > 0 & 'ligand_receptor_network' %in% names(output_list_celltype2)){
        vis_ligand_receptor_network <- output_list_celltype2[['ligand_receptor_network']]
        if(nrow(vis_ligand_receptor_network) > 0){
          plots[['ligand_receptor_network']] <- vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
        }
      }
      if(length(output_list_celltype2) > 0 & 'ligand_lfc' %in% names(output_list_celltype2)){
        vis_ligand_lfc <- output_list_celltype2[['ligand_lfc']]
        if(nrow(vis_ligand_lfc) > 0){
          plots[['ligand_lfc']] <- vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
        }
      }
      # combined the plots of the cell type
      if(length(plots) > 0){
        plot_combined <- cowplot::plot_grid(plotlist = plots)
      }
      # add resulting plot to list
      combined_plots_ct2[[cell_type2]] <- plot_combined
    }
    # add for cell type1
    combined_plots_ct1[[cell_type1]] <- combined_plots_ct2 
  }
  return(combined_plots_ct1)
}


get_lfc_celltype <- function (celltype_oi, seurat_obj, condition_colname, condition_oi, 
          condition_reference, celltype_col = "celltype", expression_pct = 0.1, logfc.threshold=0.05, test.use='wilcox') 
{
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  if (!is.null(celltype_col)) {
    seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[[celltype_col]])
    seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)
  }
  else {
    seuratObj_sender = subset(seurat_obj, idents = celltype_oi)
  }
  seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname]])
  DE_table_sender = FindMarkers(object = seuratObj_sender, 
                                ident.1 = condition_oi, ident.2 = condition_reference, 
                                min.pct = expression_pct, logfc.threshold = logfc.threshold, test.use = test.use) %>% 
    rownames_to_column("gene")
  SeuratV4 = c("avg_log2FC") %in% colnames(DE_table_sender)
  if (SeuratV4 == TRUE) {
    DE_table_sender = DE_table_sender %>% as_tibble() %>% 
      select(-p_val) %>% select(gene, avg_log2FC)
  }
  else {
    DE_table_sender = DE_table_sender %>% as_tibble() %>% 
      select(-p_val) %>% select(gene, avg_logFC)
  }
  colnames(DE_table_sender) = c("gene", celltype_oi)
  return(DE_table_sender)
}

create_dotplot_per_ct_and_tp <- function(seurat_object, features_per_ct, condition_1, condition_2, output_loc, cell_type_column='cell_type_lowerres', condition_column='timepoint.final', split_condition=F){
  Idents(seurat_object) <- cell_type_column
  # subset the Seurat object to these conditions
  seurat_object_tps <- seurat_object[, as.character(seurat_object@meta.data[[condition_column]]) == condition_1 | as.character(seurat_object@meta.data[[condition_column]]) == condition_2]
  if(split_condition){
    seurat_object_tps@meta.data$cell_type_condition <- paste(seurat_object_tps@meta.data[[condition_column]], seurat_object_tps@meta.data[[cell_type_column]], sep = '.')
    Idents(seurat_object_tps) <- 'cell_type_condition'
  }
  # check each cell type
  for(cell_type in names(features_per_ct)){
    try({
      # get the features
      features <- features_per_ct[[as.character(cell_type)]][['best_upstream_ligands']]
      print(head(features))
      print(seurat_object_tps)
      # subset the Seurat object to this cell type
      #seurat_object_tps_ct <- seurat_object_tps[, as.character(seurat_object_tps@meta.data[[cell_type_column]]) == cell_type]
      # make the plot
      p <- DotPlot(seurat_object_tps, features = features %>% rev(), cols = "RdYlBu") + RotatedAxis() + ggtitle(paste(cell_type, ' ', condition_1, ' vs ', condition_2, sep = ''))
      p
      # paste the output location together
      output_loc_full <- paste(output_loc, 'nichent_dotplot_', condition_1, '_vs_', condition_2, '_', cell_type, '.pdf', sep = '')
      ggsave(output_loc_full, width = 10, height = 10)
    })
  }
}


#
# main code
#

# where we are working
read_partition <- 'tmp01'
write_partition <- 'tmp01'

# locations of files
ligand_target_matrix_loc <- url('https://zenodo.org/record/3260758/files/ligand_target_matrix.rds')
ligand_target_matrix_loc <-paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/references/', 'ligand_target_matrix.rds', sep = '')
lr_network_loc <- url('https://zenodo.org/record/3260758/files/lr_network.rds')
lr_network_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/references/','lr_network.rds', sep = '')
weighted_networks_loc <- url("https://zenodo.org/record/3260758/files/weighted_networks.rds")
weighted_networks_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/nichenet/references/','weighted_networks.rds', sep = '')
objects_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/', sep = '')
combined_v2_loc <- paste(objects_loc, 'combined.v2.20210629.ct.rds', sep = '')
combined_v3_loc <- paste(objects_loc, 'combined.v3.20210629.ct.rds', sep = '')

# read network data
ligand_target_matrix = readRDS(ligand_target_matrix_loc)
lr_network = readRDS(lr_network_loc)
weighted_networks = readRDS(weighted_networks_loc)

# read the object
combined_v2 <- readRDS(combined_v2_loc)
combined_v2 <- combined_v2[, combined_v2@meta.data$cell_type_lowerres %in% c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')]
# compare Baseline to t24h
v2_Baseline_vs_t24h <- do_nichenet_analysis_per_celltype(combined_v2, 'timepoint.final', 't24h', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v2_Baseline_vs_t24h, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v2_Baseline_vs_t24h_nichenet_onlymajors.rds')
# compare Baseline to t8w
v2_Baseline_vs_t8w <- do_nichenet_analysis_per_celltype(combined_v2, 'timepoint.final', 't8w', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v2_Baseline_vs_t8w, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v2_Baseline_vs_t8w_nichenet_onlymajor.rds')
# compare HC to Baseline
v2_UT_vs_Baseline <- do_nichenet_analysis_per_celltype(combined_v2, 'timepoint.final', 'Baseline', 'UT', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v2_UT_vs_Baseline, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v2_UT_vs_Baseline_nichenet_onlymajor.rds')

# trying to create the plots
v2_Baseline_vs_t24h_plots <- nichenet_output_to_plot(v2_Baseline_vs_t24h)
v2_Baseline_vs_t8w_plots <- nichenet_output_to_plot(v2_Baseline_vs_t8w)
v2_UT_vs_Baseline_plots <- nichenet_output_to_plot(v2_UT_vs_Baseline)

# save plots
for(plot_name in names(v2_Baseline_vs_t24h_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined_onlymajor/', 'nichenet_', 'Baseline_vs_t24h_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v2_Baseline_vs_t8w_plots[[plot_name]])
}
for(plot_name in names(v2_Baseline_vs_t8w_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined_onlymajor/', 'nichenet_', 'Baseline_vs_t8w_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v2_Baseline_vs_t8w_plots[[plot_name]])
}
for(plot_name in names(v2_UT_vs_Baseline_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined_onlymajor/', 'nichenet_', 'UT_vs_Baseline_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v2_UT_vs_Baseline_plots[[plot_name]])
}


# read the object
combined_v3 <- readRDS(combined_v3_loc)
combined_v3 <- combined_v3[, combined_v3@meta.data$cell_type_lowerres %in% c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')]
# compare Baseline to t24h
v3_Baseline_vs_t24h <- do_nichenet_analysis_per_celltype(combined_v3, 'timepoint.final', 't24h', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v3_Baseline_vs_t24h, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v3_Baseline_vs_t24h_nichenet_onlymajor.rds')
# compare Baseline to t8w
v3_Baseline_vs_t8w <- do_nichenet_analysis_per_celltype(combined_v3, 'timepoint.final', 't8w', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v3_Baseline_vs_t8w, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v3_Baseline_vs_t8w_nichenet_onlymajor.rds')
# compare HC to Baseline
v3_UT_vs_Baseline <- do_nichenet_analysis_per_celltype(combined_v3, 'timepoint.final', 'Baseline', 'UT', lr_network, weighted_networks, ligand_target_matrix)
saveRDS(v3_UT_vs_Baseline, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/v3_UT_vs_Baseline_nichenet_onlymajor.rds')

# trying to create the plots
v3_Baseline_vs_t24h_plots <- nichenet_output_to_plot(v3_Baseline_vs_t24h)
v3_Baseline_vs_t8w_plots <- nichenet_output_to_plot(v3_Baseline_vs_t8w)
v3_UT_vs_Baseline_plots <- nichenet_output_to_plot(v3_UT_vs_Baseline)

# save plots
for(plot_name in names(v3_Baseline_vs_t24h_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor//', 'nichenet_', 'Baseline_vs_t24h_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v3_Baseline_vs_t24h_plots[[plot_name]])
}
for(plot_name in names(v3_Baseline_vs_t8w_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor//', 'nichenet_', 'Baseline_vs_t8w_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v3_Baseline_vs_t8w_plots[[plot_name]])
}
for(plot_name in names(v3_UT_vs_Baseline_plots)){
  ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor//', 'nichenet_', 'UT_vs_Baseline_', plot_name, '_onlymajor.pdf', sep = ''), width = 20, heigh = 20, plot=v3_UT_vs_Baseline_plots[[plot_name]])
}

# do specifically for one cell type now
v2_Baseline_vs_t24h_perct <- do_nichenet_analysis_versus_each_celltype(combined_v2, 'timepoint.final', 't24h', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
v2_Baseline_vs_t8w_perct <- do_nichenet_analysis_versus_each_celltype(combined_v2, 'timepoint.final', 't8w', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
v2_t8w_vs_Baseline_perct <- do_nichenet_analysis_versus_each_celltype(combined_v2, 'timepoint.final', 'Baseline', 'UT', lr_network, weighted_networks, ligand_target_matrix)

# save the plots once again
for(ct1 in names(v2_Baseline_vs_t24h_perct)){
  for(ct2 in names(v2_Baseline_vs_t24h_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t24h_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v2_Baseline_vs_t24h_perct[[ct1]][[ct2]])
  }
}
# save the plots once again
for(ct1 in names(v2_Baseline_vs_t8wh_perct)){
  for(ct2 in names(v2_Baseline_vs_t8w_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t8w_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v2_Baseline_vs_t8w_perct[[ct1]][[ct2]])
  }
}
# save the plots once again
for(ct1 in names(v2_Baseline_vs_t8wh_perct)){
  for(ct2 in names(v2_Baseline_vs_t8w_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t8w_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v2_Baseline_vs_t8w_perct[[ct1]][[ct2]])
  }
}

# do specifically for one cell type now
v3_Baseline_vs_t24h_perct <- do_nichenet_analysis_versus_each_celltype(combined_v3, 'timepoint.final', 't24h', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
v3_Baseline_vs_t8w_perct <- do_nichenet_analysis_versus_each_celltype(combined_v3, 'timepoint.final', 't8w', 'Baseline', lr_network, weighted_networks, ligand_target_matrix)
v3_t8w_vs_Baseline_perct <- do_nichenet_analysis_versus_each_celltype(combined_v3, 'timepoint.final', 'Baseline', 'UT', lr_network, weighted_networks, ligand_target_matrix)

# save the plots once again
for(ct1 in names(v3_Baseline_vs_t24h_perct)){
  for(ct2 in names(v3_Baseline_vs_t24h_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t24h_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v3_Baseline_vs_t24h_perct[[ct1]][[ct2]])
  }
}
# save the plots once again
for(ct1 in names(v3_Baseline_vs_t8wh_perct)){
  for(ct2 in names(v3_Baseline_vs_t8w_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t8w_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v3_Baseline_vs_t8w_perct[[ct1]][[ct2]])
  }
}
# save the plots once again
for(ct1 in names(v3_Baseline_vs_t8wh_perct)){
  for(ct2 in names(v3_Baseline_vs_t8w_perct[[ct1]])){
    ggsave(paste('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_per_ct/', 'nichenet_', 'Baseline_vs_t8w_', ct1, '_vs_', ct2, '.pdf', sep = ''), width = 20, heigh = 20, plot=v3_Baseline_vs_t8w_perct[[ct1]][[ct2]])
  }
}




create_dotplot_per_ct_and_tp(combined_v3, v3_Baseline_vs_t24h, 'Baseline', 't24h', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor/', split_condition = T)
create_dotplot_per_ct_and_tp(combined_v3, v3_Baseline_vs_t8w, 'Baseline', 't8w', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor/', split_condition = T)
create_dotplot_per_ct_and_tp(combined_v3, v3_UT_vs_Baseline, 'UT', 'Baseline', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v3_combined_onlymajor/', split_condition = T)
create_dotplot_per_ct_and_tp(combined_v2, v2_Baseline_vs_t24h, 'Baseline', 't24h', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined_onlymajor/', split_condition = T)
create_dotplot_per_ct_and_tp(combined_v2, v2_Baseline_vs_t8w, 'Baseline', 't8w', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined_onlymajor/', split_condition = T)
create_dotplot_per_ct_and_tp(combined_v2, v2_UT_vs_Baseline, 'UT', 'Baseline', '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/nichenet/MAST/v2_combined/_onlymajor', split_condition = T)


