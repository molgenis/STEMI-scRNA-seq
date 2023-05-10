#
# description of script
#

#################
# libraries     #
#################

library(Seurat)
library(Matrix)
library(VariantAnnotation)
library(ggplot2)
library(UCell)

#################
# functions     #
#################

add_pathway_module_score <- function(seurat_object, pathway_gene_list_loc, assay='SCT', name='pathway'){
  # read the file, only one column so we can grab that straight away
  pathway_genes <- read.table(pathway_gene_list_loc, header = F)$V1
  # add the module score
  seurat_object <- AddModuleScore(seurat_object, list(pathway_genes), assay = assay, name = name)
  # seurat appends a '1' to the module score, fix that
  seurat_object@meta.data[[name]] <- seurat_object@meta.data[[paste(name, '1', sep = '')]]
  seurat_object@meta.data[[paste(name, '1', sep = '')]] <- NULL
  return(seurat_object)
}

add_pathway_ucell_score <- function(seurat_object, pathway_gene_list_loc, assay='SCT', name='pathway'){
  # read the file, only one column so we can grab that straight away
  pathway_genes <- read.table(pathway_gene_list_loc, header = F)$V1
  # add the module score
  pathway_genes_list <- list()
  pathway_genes_list[[name]] <- pathway_genes
  seurat_object <- AddModuleScore_UCell(seurat_object, pathway_genes_list, assay = assay, name = '')
  # seurat appends a '_UCell' to the module score, fix that
  #seurat_object@meta.data[[name]] <- seurat_object@meta.data[[paste(name, '_UCell', sep = '')]]
  #seurat_object@meta.data[[paste(name, '_UCell', sep = '')]] <- NULL
  return(seurat_object)
}


add_pathway_scores <- function(seurat_object, pathway_gene_locations_per_pathway, assay='SCT', use_ucell=F){
  # check each named location
  for(pathway_name in names(pathway_gene_locations_per_pathway)){
    # grab the specific pathway location
    pathway_loc <- pathway_gene_locations_per_pathway[[pathway_name]]
    # do the scoring
    if(use_ucell){
      seurat_object <- add_pathway_ucell_score(seurat_object, pathway_loc, assay = assay, name = pathway_name)
    }
    else{
      seurat_object <- add_pathway_module_score(seurat_object, pathway_loc, assay = assay, name = pathway_name)
    }
  }
  return(seurat_object)
}

get_mean_expression_matrices_metadata <- function(seurat_object, metadata_columns, cell_type_column='cell_type', condition_column='timepoint.final', assignment_column='assignment.final'){
  # we need to store the output somewhere
  mean_expression_per_condition <- list()
  # check each condition
  for(condition in as.character(unique(seurat_object@meta.data[[condition_column]]))){
    print(condition)
    # we need to store per cell type as well
    mean_expression_per_condition_cell_type <- list()
    # subset to that condition
    seurat_object_condition <- seurat_object[, as.character(seurat_object@meta.data[[condition_column]]) == condition]
    # check each cell type
    for(cell_type in as.character(unique(seurat_object_condition@meta.data[[cell_type_column]]))){
      print(cell_type)
      # subset to that cell type
      seurat_object_condition_celltype <- seurat_object_condition[, as.character(seurat_object_condition@meta.data[[cell_type_column]]) == cell_type]
      # get the unique participants in this subset
      participants <- as.character(unique(seurat_object_condition_celltype@meta.data[[assignment_column]]))
      # remove NAs
      participants <- participants[!is.na(participants)]
      # create an expression matrix
      expression_matrix <- matrix(, ncol = length(participants), nrow = length(metadata_columns), dimnames = list(metadata_columns, participants))
      # check each participant
      for(participant in participants){
        # check each metadata column
        for(metadata_column in metadata_columns){
          # calculate the mean
          mean_value_metadata_column <- mean(seurat_object_condition_celltype@meta.data[as.character(seurat_object_condition_celltype@meta.data[[assignment_column]]) == participant, metadata_column])
          # add to the matrix
          expression_matrix[metadata_column, participant] <- mean_value_metadata_column
        }
      }
      # add info of this cell type
      mean_expression_per_condition_cell_type[[cell_type]] <- expression_matrix
    }
    # add to data per condition
    mean_expression_per_condition[[condition]] <- mean_expression_per_condition_cell_type
  }
  return(mean_expression_per_condition)
}

plottable_pathways_from_metadata <- function(seurat_object, pathway_columns, mean_first=T, cell_type_column='cell_type', condition_column='timepoint.final', assignment_column='assignment.final'){
  full_plot_data <- NULL
  # the mean method means the score per participant
  if(mean_first){
    # check each pathway
    table_per_condition_celltype <- get_mean_expression_matrices_metadata(seurat_object, pathway_columns, cell_type_column = cell_type_column, condition_column = condition_column, assignment_column = assignment_column)
    # go through each condition
    for(condition in names(table_per_condition_celltype)){
      # go through each cell type as well
      for(cell_type in names(table_per_condition_celltype[[condition]])){
        # get table of this cell type and condition
        table_cell_and_condition <- table_per_condition_celltype[[condition]][[cell_type]]
        # get the participants
        participants <- colnames(table_cell_and_condition)
        # transpose, so that the pathways are columns, and participants rows
        table_cell_and_condition <- data.frame(t(table_cell_and_condition))
        # add participant as column as well
        table_cell_and_condition$participant <- participants
        # add the condition and cell type data
        table_cell_and_condition$cell_type <- cell_type
        table_cell_and_condition$condition <- condition
        # check if this is the first round
        if(is.null(full_plot_data)){
          full_plot_data <- table_cell_and_condition
        }
        else{
          full_plot_data <- rbind(full_plot_data, table_cell_and_condition)
        }
      }
    }
  }
  else{
    # get the data
    full_plot_data <- seurat_object@meta.data[, c(pathway_columns, assignment_column, cell_type_column, condition_column)]
    # rename some of the column
    if(assignment_column != 'participant'){
      full_plot_data$participant <- full_plot_data[[assignment_column]]
      full_plot_data[[assignment_column]] <- NULL
    }
    if(cell_type_column != 'cell_type'){
      full_plot_data$cell_type <- full_plot_data[[cell_type_column]]
      full_plot_data[[cell_type_column]] <- NULL
    }
    if(condition_column != 'condition'){
      full_plot_data$condition <- full_plot_data[[condition_column]]
      full_plot_data[[condition_column]] <- NULL
    }
  }
  return(full_plot_data)
}


plot_pathways <- function(plottable_pathways, pathways_to_plot, output_loc, use_label_dict = T, use_color_dict=T, pointless=F, legendless=F, angle_x_labels=F, paper_style=F, width=10, height=10){
  # if requested, clean up the labels
  if(use_label_dict){
    #plottable_pathways$cell_type_nicer <- as.vector(unlist(label_dict()[as.character(plottable_pathways$cell_type)]))
    plottable_pathways$condition <- as.vector(unlist(label_dict()[as.character(plottable_pathways$condition)]))
  }
  # check each cell type
  for(cell_type in unique(plottable_pathways[['cell_type']])){
    # subset to cell type
    plottable_cell_type <- plottable_pathways[plottable_pathways[['cell_type']] == cell_type, ]
    # make plot for each pathway
    for(pathway in intersect(pathways_to_plot, colnames(plottable_cell_type))){
      # set as a column
      plottable_cell_type$score <- plottable_cell_type[[pathway]]
      # create the plot
      p <- ggplot(data = plottable_cell_type, mapping = aes(x = condition, y = score, fill=condition)) + geom_boxplot(outlier.shape = NA)
      # use custom colours if requesteds
      if(use_color_dict){
        # grab colours for each condition
        colFill <- scale_fill_manual(name = 'condition', values = unlist(get_color_coding_dict()[as.character(plottable_cell_type$condition)]))
        # add to plot
        p <- p + colFill
      }
      # other options
      if(pointless){
        p <- p + theme(axis.text.x=element_blank(), 
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank())
      }
      if(legendless){
        p <- p + theme(legend.position = 'none')
      }
      if(paper_style){
        p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
      }
      if(angle_x_labels){
        p <- p + theme(axis.text.x = element_text(angle = 90))
      }
      if(use_label_dict){
        p <- p + ggtitle(paste(pathway, 'score', 'in', label_dict()[[cell_type]]))
      }
      else{
        p <- p + ggtitle(paste(pathway, 'score', 'in', cell_type))
      }
      # paste together the output location
      output_loc_full <- paste(output_loc, 'pathways_', pathway, '_', cell_type, '.pdf', sep = '')
      print(output_loc_full)
      ggsave(filename = output_loc_full, plot = p, width = width, height = height)
    }
  }
}


plot_qtl <- function(genotype_data, expression_data, snp, y, xlab='genotype', ylab='y', paper_style=F){
  # get the overlapping samples
  samples_to_use <- intersect(colnames(genotype_data), colnames(expression_data))
  # grab the data from the two matrices
  genotypes <- genotype_data[snp, samples_to_use]
  y_values <- expression_data[y, samples_to_use]
  # turn into a plottable frame
  plot_data <- data.frame(genotype=genotypes, y=y_values)
  # do the plotting
  p <- ggplot(data = plot_data, mapping = aes(x = genotype, y=y, fill=genotype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.2) +
    ggtitle('')
  if(length(unique(genotypes)) == 1){
    p <- p + scale_fill_manual(values=c("#57a350"), guide = F)
  }
  else if(length(unique(genotypes)) == 2){
    p <- p + scale_fill_manual(values=c("#57a350", "#fd7600"), guide = F)
  }
  else if(length(unique(genotypes)) == 3){
    print('manual colors')
    p <- p + scale_fill_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F)
  }
  if(paper_style){
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["t0-HC"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["t24h-HC"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t8w-HC"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t24h-t0"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t8w-t0"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
  color_coding[["t8w-t24h"]] <- "#C00040" #"rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  # this one should not exist
  color_coding[['doublet']] <- #4287f5
  color_coding[['unclassified']] <- '#dbdbdb'
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[['UTBaseline']] <- 'UT-Baseline'
  label_dict[['UTt24h']] <- 'UT-t24h'
  label_dict[['UTt8w']] <- 'UT-t8w'
  label_dict[['Baselinet24h']] <- 'Baseline-t24h'
  label_dict[['Baselinet8w']] <- 'Baseline-t8w'
  label_dict[['t24ht8w']] <- 't24h-t8w'
  label_dict[['UTBaseline']] <- 'HC-t0'
  label_dict[['UTt24h']] <- 'HC-t24h'
  label_dict[['UTt8w']] <- 'HC-t8w'
  label_dict[['Baselinet24h']] <- 't0-t24h'
  label_dict[['Baselinet8w']] <- 't0-t8w'
  # conditions
  label_dict[['UT']] <- 'HC'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  label_dict[['doublet']] <- 'doublet'
  label_dict[['unclassified']] <- 'unclassified'
  return(label_dict)
}


#################
# main code     #
#################

# partition for reading
read_partition <- 'tmp01/'
write_partition <- 'tmp01/'

# locations of files
pathway_annotation_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/annotations/', sep = '')
pathway_annotation_il1 <- paste(pathway_annotation_loc, 'il1-pathway-reactome.txt', sep = '')
pathway_annotation_il6 <- paste(pathway_annotation_loc, 'il6-pathway-reactome.txt', sep = '')
pathway_annotation_il10 <- paste(pathway_annotation_loc, 'il10-pathway-reactome.txt', sep = '')
pathway_annotation_il17 <- paste(pathway_annotation_loc, 'il17-pathway-reactome.txt', sep = '')
pathway_annotation_il23 <- paste(pathway_annotation_loc, 'il23-pathway-reactome.txt', sep = '')
pathway_annotation_il4il13 <- paste(pathway_annotation_loc, 'il4il13-pathway-reactome.txt', sep = '')
# put into a list
pathways_list <- list(il1=pathway_annotation_il1, il6=pathway_annotation_il6, il10=pathway_annotation_il6, il17=pathway_annotation_il17, il23=pathway_annotation_il23, il4il13=pathway_annotation_il4il13)

# location of genotype data
genotypes_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi.vcf.gz', sep = '')

# location of plots
pathway_plots_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/plots/pathway_scores/', sep = '')
module_score_pathway_plots_loc <- paste(pathway_plots_loc, 'seurat_module_score/', sep = '')
v2_module_score_pathway_plots_loc <- paste(module_score_pathway_plots_loc, 'v2/', sep = '')
v3_module_score_pathway_plots_loc <- paste(module_score_pathway_plots_loc, 'v3/', sep = '')
ucell_score_pathway_plots_loc <- paste(pathway_plots_loc, 'ucell_score/', sep = '')
v2_ucell_score_pathway_plots_loc <- paste(ucell_score_pathway_plots_loc, 'v2/', sep = '')
v3_ucell_score_pathway_plots_loc <- paste(ucell_score_pathway_plots_loc, 'v3/', sep = '')

# location of Seurat things
object_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/', sep = '')
cardio.integrated_loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')
v2_combined_loc <- paste(object_loc, 'combined.v2.20210629.ct.rds', sep = '')
v3_combined_loc <- paste(object_loc, 'combined.v3.20210629.ct.rds', sep = '')

# load objects
v2_combined <- readRDS(v2_combined_loc)
v3_combined <- readRDS(v3_combined_loc)

# add pathways
v2_combined <- add_pathway_scores(v2_combined, pathways_list)
v3_combined <- add_pathway_scores(v3_combined, pathways_list)
# get plottable frame
v2_combined_plottable <- plottable_pathways_from_metadata(v2_combined, names(pathways_list), cell_type_column = 'cell_type_lowerres')
v3_combined_plottable <- plottable_pathways_from_metadata(v3_combined, names(pathways_list), cell_type_column = 'cell_type_lowerres')
# plot each celltype per condition
plot_pathways(v2_combined_plottable, names(pathways_list), v2_module_score_pathway_plots_loc)
plot_pathways(v3_combined_plottable, names(pathways_list), v3_module_score_pathway_plots_loc)

# do the same thing for the ucell ones
v2_combined <- add_pathway_scores(v2_combined, pathways_list, use_ucell = T)
v3_combined <- add_pathway_scores(v3_combined, pathways_list, use_ucell = T)
v2_combined_plottable <- plottable_pathways_from_metadata(v2_combined, names(pathways_list), cell_type_column = 'cell_type_lowerres')
v3_combined_plottable <- plottable_pathways_from_metadata(v3_combined, names(pathways_list), cell_type_column = 'cell_type_lowerres')
plot_pathways(v2_combined_plottable, names(pathways_list), v2_ucell_score_pathway_plots_loc)
plot_pathways(v3_combined_plottable, names(pathways_list), v3_ucell_score_pathway_plots_loc)


v2_combined_il6_features <- get_mean_expression_matrices_metadata(v2_combined, c('IL6'), cell_type_column = 'cell_type_lowerres')

# load genotype data
genotypes <- readGT(genotypes_loc)
genotypes['rs6689306', genotypes['rs6689306', ] == '1/0'] <- '0/1'
genotypes['rs6689306', genotypes['rs6689306', ] == '1|0'] <- '0/1'
genotypes['rs6689306', genotypes['rs6689306', ] == '0|1'] <- '0/1'
genotypes['rs6689306', genotypes['rs6689306', ] == '1|1'] <- '1/1'
genotypes['rs6689306', genotypes['rs6689306', ] == '0|0'] <- '0/0'

# load genotype data
genotypes <- readGT(genotypes_loc)
genotypes['rs7529229', genotypes['rs7529229', ] == '1/0'] <- '0/1'
genotypes['rs7529229', genotypes['rs7529229', ] == '1|0'] <- '0/1'
genotypes['rs7529229', genotypes['rs7529229', ] == '0|1'] <- '0/1'
genotypes['rs7529229', genotypes['rs7529229', ] == '1|1'] <- '1/1'
genotypes['rs7529229', genotypes['rs7529229', ] == '0|0'] <- '0/0'

stemi_CD4T_Baseline_rs7529229_IL6_Ucell <- plot_qtl(genotype=genotypes, v2_combined_il6_ucell_features[['Baseline']][['CD4T']], 'rs7529229', 'IL6_Ucell', 'pathway score', T) + ggtitle('v2 Baseline')
stemi_CD4T_UT_rs7529229_IL6_Ucell <- plot_qtl(genotype=genotypes, v2_combined_il6_ucell_features[['UT']][['CD4T']], 'rs7529229', 'IL6_Ucell', 'pathway score', T) + ggtitle('v2 UT')
stemi_CD4T_t24h_rs7529229_IL6_Ucell <- plot_qtl(genotype=genotypes, v2_combined_il6_ucell_features[['t24h']][['CD4T']], 'rs7529229', 'IL6_Ucell', 'pathway score', T) + ggtitle('v2 t24h')
stemi_CD4T_t8w_rs7529229_IL6_Ucell <- plot_qtl(genotype=genotypes, v2_combined_il6_ucell_features[['t8w']][['CD4T']], 'rs7529229', 'IL6_Ucell', 'pathway score', T) + ggtitle('v2 t8w')
plot_grid(plotlist = list(stemi_CD4T_UT_rs7529229_IL6_Ucell, stemi_CD4T_Baseline_rs7529229_IL6_Ucell, stemi_CD4T_t24h_rs7529229_IL6_Ucell, stemi_CD4T_t8w_rs7529229_IL6_Ucell))
ggsave('stemi_CD4T_v2_rs7529229_IL6_Ucell_pathway.pdf', width = 10, height = 10)



