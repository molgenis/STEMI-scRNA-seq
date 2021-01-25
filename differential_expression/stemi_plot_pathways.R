####################
# libraries        #
####################

library(Seurat)
library(ggplot2)


####################
# Functions        #
####################

# plot the module score in a boxplot
plot_average_expression <- function(seurat_object, module_score_column_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', title='pathway', color_by_ct=T){
  # get the metadata
  metadata <- seurat_object@meta.data
  # limit the cell types
  metadata <- metadata[metadata[[cell_type_column]] %in% cell_types, ]
  # limit the conditions
  metadata <- metadata[metadata[[condition_column]] %in% conditions, ]
  # set some hardcoded column names
  metadata$ct_column <- metadata[[cell_type_column]]
  metadata$tp_column <- metadata[[condition_column]]
  metadata$msc_column <- metadata[[module_score_column_name]]
  print(head(metadata))
  metadata$cttp_column <- paste(metadata$ct_column, metadata$tp_column, sep = '\n')
  # grab our dict for colors
  cc <- get_color_coding_dict()
  # if the colouring is by cell type
  if(color_by_ct){
    colScale <- scale_fill_manual(name = metadata$ct_column, values = unlist(cc[cell_types]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=ct_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition')
  }
  # if the colouring is by the condition
  else{
    ct_tp_order <- c()
    for(cell_type in cell_types){
      ct_tp_order <- c(ct_tp_order, paste(rep(cell_type, each = length(conditions)), conditions, sep = "\n"))
    }
    colScale <- scale_fill_manual(name = metadata$tp_column, values = unlist(cc[conditions]))
    ggplot(metadata, aes(x=cttp_column, y=msc_column, fill=tp_column)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      labs(y = 'module scores', x='cell type and condition') +
      scale_x_discrete(limits = ct_tp_order)
  }
}

# add the module score from a table of genes (now using as pathway plot)
add_module_score_from_table <- function(pathway_gene_table_loc, pathway_name, seurat_object, and_plot=T, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = T){
  # get the cytokine genes
  pathway_df <- read.table(pathway_gene_table_loc)
  pathway_genes <- pathway_df$V1
  pathway_genes_in_seurat_object <- intersect(rownames(seurat_object), pathway_genes)
  # add to list
  pathway_list_seurat_object <- list()
  pathway_list_seurat_object[[pathway_name]] <- pathway_genes_in_seurat_object
  # add the module score
  seurat_object <- AddModuleScore(seurat_object, features = pathway_list_seurat_object, name = pathway_name)
  if(and_plot){
    plot_average_expression(seurat_object, module_score_column_name=paste(pathway_name, '1', sep = ''), cell_types=cell_types, conditions=conditions, cell_type_column=cell_type_column, condition_column=condition_column, title=pathway_name, color_by_ct = color_by_ct)
  }
  return(seurat_object)
}

# mapping of conditions and cell types to colours
get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- 'grey'
  color_coding[["3hCA"]] <- "khaki2"
  color_coding[["24hCA"]] <- "khaki4"
  color_coding[["3hMTB"]] <- "paleturquoise1"
  color_coding[["24hMTB"]] <- "paleturquoise3"
  color_coding[["3hPA"]] <- "rosybrown1"
  color_coding[["24hPA"]] <- "rosybrown3"
  color_coding[["X3hCA"]] <- "khaki2"
  color_coding[["X24hCA"]] <- "khaki4"
  color_coding[["X3hMTB"]] <- "paleturquoise1"
  color_coding[["X24hMTB"]] <- "paleturquoise3"
  color_coding[["X3hPA"]] <- "rosybrown1"
  color_coding[["X24hPA"]] <- "rosybrown3"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}

####################
# Main code        #
####################

# set location of object
cardio.integrated.loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated_20201209.rds'
cardio.integrated <- readRDS(cardio.integrated.loc)

# location of the pathway genes
pathway_gene_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/annotations/'
pathway_files <- c('KEGG_IL-17_signaling_pathway.txt', 'KEGG_Th17_cell_differentiation.txt', 'MySigDBC2BIOCARTA_IL23-mediated_signaling_events.txt', 'MySigDBC2BIOCARTA_IL6-mediated_signaling_events.txt', 'REACTOME_Interleukin-10_signaling.txt', 'REACTOME_Interleukin-4_and_13_signalling.txt')

# location to output the plots
plot_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/plots/'

# subset to v2 and v3
cardio.v2 <- subset(cardio.integrated, subset = chem == 'V2')
cardio.v3 <- subset(cardio.integrated, subset = chem == 'V3')

# plot each pathway
for(pathway_file in pathway_files){
  # set the full file loc of the gene list
  full_pathway_file <- paste(pathway_gene_loc, pathway_file, sep='')
  # create a regex to get the last index of .
  last_dot_pos <- "\\.[^\\.]*$"
  # this allows us to remove the filename extention, which we will use as the name
  pathway_name <- substr(file, 1, regexpr(last_dot_pos, pathway_file)-1)
  # do for v2
  add_module_score_from_table(full_pathway_file, pathway_name, cardio.v2, cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F)
  output_loc_v2 <- paste(plot_output_loc, 'mono_v2_', pathway_name, '.png', sep = '')
  ggsave(output_loc_v2, width = 10, height = 10)
  # and for v3
  add_module_score_from_table(full_pathway_file, pathway_name, cardio.v3, cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F)
  output_loc_v2 <- paste(plot_output_loc, 'mono_v3_', pathway_name, '.png', sep = '')
  ggsave(output_loc_v2, width = 10, height = 10)
}


