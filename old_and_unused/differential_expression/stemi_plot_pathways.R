####################
# libraries        #
####################

library(Seurat)
library(ggplot2)


####################
# Functions        #
####################

# plot the module score in a boxplot
plot_average_expression <- function(seurat_object, module_score_column_name, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', title='pathway', color_by_ct=T, to_per_part=F, participant_column='assignment.final'){
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
  metadata$cttp_column <- paste(metadata$ct_column, metadata$tp_column, sep = '\n')
  # convert to per_participant if requested
  if(to_per_part){
    metadata <- metadata_to_average_module_per_part(metadata, cttp_column='cttp_column', msc_column='msc_column', ct_column='ct_column', condition_column='tp_column', participant_column=participant_column)
  }
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

metadata_to_average_module_per_part <- function(metadata, cttp_column='cttp_column', msc_column='msc_column', ct_column='ct_column', condition_column='tp_column', participant_column='assignment.final'){
  metadata_averaged <- NULL
  for(cttp in unique(metadata[[cttp_column]])){
    # subset to cell type and condition combo
    metadata_cttp <- metadata[metadata[[cttp_column]] == cttp ,]
    # check for each participant here
    for(participant in unique(metadata_cttp[[participant_column]])){
      # calculate the mean
      mean_participant <- mean(metadata_cttp[metadata_cttp[[participant_column]] == participant, msc_column])
      # get other data
      ct <- as.character(metadata_cttp[metadata_cttp[[participant_column]] == participant, ct_column][1])
      tp <- as.character(metadata_cttp[metadata_cttp[[participant_column]] == participant, condition_column][1])
      # set as a dataframe
      metadata_averaged_row <- data.frame(cttp_column=c(cttp), msc_column=c(mean_participant), ct_column=c(ct), condition_column=c(tp), participant_column=c(participant), stringsAsFactors = F)
      # set as initial or append
      if(is.null(metadata_averaged)){
        metadata_averaged <- metadata_averaged_row
      }
      else{
        metadata_averaged <- rbind(metadata_averaged, metadata_averaged_row)
      }
    }
  }
  # set colnames as before
  colnames(metadata_averaged) <- c(cttp_column, msc_column, ct_column, condition_column, participant_column)
  return(metadata_averaged)
}

# add the module score from a table of genes (now using as pathway plot)
add_module_score_from_table <- function(pathway_gene_table_loc, pathway_name, seurat_object, and_plot=T, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = T, to_per_part=F, participant_column='assignment.final'){
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
    plot_average_expression(seurat_object, module_score_column_name=paste(pathway_name, '1', sep = ''), cell_types=cell_types, conditions=conditions, cell_type_column=cell_type_column, condition_column=condition_column, title=pathway_name, color_by_ct = color_by_ct, to_per_part=to_per_part, participant_column=participant_column)
  }
  return(seurat_object)
}

# mapping of conditions and cell types to colours
get_color_coding_dict <- function(){
  # set the condition combo colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  # set the condition colors
  color_coding[["UT"]] <- "#D4ECF4"
  color_coding[["t24h"]] <- "#008BC4"
  color_coding[["t8w"]] <- "#262D60"
  color_coding[["Baseline"]] <- "#66C5DD"
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
cardio.integrated.loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20201209.rds'
cardio.integrated <- readRDS(cardio.integrated.loc)

# location of the pathway genes
pathway_gene_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/annotations/'
pathway_files <- c('KEGG_IL-17_signaling_pathway.txt', 'KEGG_Th17_cell_differentiation.txt', 'MySigDBC2BIOCARTA_IL23-mediated_signaling_events.txt', 'MySigDBC2BIOCARTA_IL6-mediated_signaling_events.txt', 'REACTOME_Interleukin-10_signaling.txt', 'REACTOME_Interleukin-4_and_13_signalling.txt')

# location to output the plots
plot_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/plots/'

# subset to v2 and v3
cardio.v2 <- subset(cardio.integrated, subset = chem == 'V2')
cardio.v3 <- subset(cardio.integrated, subset = chem == 'V3')

# normalize RNA
DefaultAssay(cardio.v2) <- 'RNA'
cardio.v2 <- NormalizeData(cardio.v2)
DefaultAssay(cardio.v3) <- 'RNA'
cardio.v3 <- NormalizeData(cardio.v3)
# normalize SCT
DefaultAssay(cardio.v2) <- 'SCT'
cardio.v2 <- SCTransform(cardio.v2)
DefaultAssay(cardio.v3) <- 'SCT'
cardio.v3 <- SCTransform(cardio.v3)

# plot each pathway
for(pathway_file in pathway_files){
  # set the full file loc of the gene list
  full_pathway_file <- paste(pathway_gene_loc, pathway_file, sep='')
  # create a regex to get the last index of .
  last_dot_pos <- "\\.[^\\.]*$"
  # this allows us to remove the filename extention, which we will use as the name
  pathway_name <- substr(pathway_file, 1, regexpr(last_dot_pos, pathway_file)-1)
  # we need to replace dashes with underscores
  pathway_name <- gsub('-', '_', pathway_name)
  # do for each cell type
  for(cell_type in unique(cardio.integrated@meta.data$cell_type_lowerres)){
    # do for v2
    DefaultAssay(cardio.v2) <- 'RNA'
    add_module_score_from_table(full_pathway_file, pathway_name, cardio.v2, cell_types=c(cell_type), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F, to_per_part=T, participant_column='assignment.final')
    output_loc_v2 <- paste(plot_output_loc, cell_type, '_v2_', pathway_name, '_RNA_perpart.png', sep = '')
    ggsave(output_loc_v2, width = 10, height = 10)
    # do for v2
    DefaultAssay(cardio.v2) <- 'SCT'
    add_module_score_from_table(full_pathway_file, pathway_name, cardio.v2, cell_types=c(cell_type), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F, to_per_part=T, participant_column='assignment.final')
    output_loc_v2 <- paste(plot_output_loc, cell_type, '_v2_', pathway_name, '_SCT_perpart.png', sep = '')
    ggsave(output_loc_v2, width = 10, height = 10)
    # and for v3
    DefaultAssay(cardio.v3) <- 'RNA'
    add_module_score_from_table(full_pathway_file, pathway_name, cardio.v3, cell_types=c(cell_type), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F, to_per_part=T, participant_column='assignment.final')
    output_loc_v2 <- paste(plot_output_loc, cell_type, '_v3_', pathway_name, '_RNA_perpart.png', sep = '')
    ggsave(output_loc_v2, width = 10, height = 10)
    # and for v3
    DefaultAssay(cardio.v3) <- 'SCT'
    add_module_score_from_table(full_pathway_file, pathway_name, cardio.v3, cell_types=c(cell_type), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = F, to_per_part=T, participant_column='assignment.final')
    output_loc_v2 <- paste(plot_output_loc, cell_type, '_v3_', pathway_name, '_SCT_perpart.png', sep = '')
    ggsave(output_loc_v2, width = 10, height = 10)
  }
}


