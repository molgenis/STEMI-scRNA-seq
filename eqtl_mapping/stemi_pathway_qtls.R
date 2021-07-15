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
  seurat_object <- AddModuleScore_UCell(seurat_object, list(pathway = pathway_genes), assay = assay, name = name)
  # seurat appends a '_UCell' to the module score, fix that
  seurat_object@meta.data[[name]] <- seurat_object@meta.data[[paste(name, '_UCell', sep = '')]]
  seurat_object@meta.data[[paste(name, '_UCell', sep = '')]] <- NULL
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


#################
# main code     #
#################

# partition for reading
read_partition <- 'tmp04/'
write_partition <- 'tmp04/'

# locations of files
pathway_annotation_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/pathways/annotations/', sep = '')
pathway_annotation_il6 <- paste(pathway_annotation_loc, 'il6-pathway-reactome.txt', sep = '')

# location of genotype data
genotypes_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi.vcf.gz', sep = '')

# location of Seurat things
object_loc <- paste('/groups/umcg-wijmenga/', read_partition, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/', sep = '')
cardio.integrated_loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')
v2_combined_loc <- paste(object_loc, 'combined.v2.20210629.ct.rds', sep = '')
v3_combined_loc <- paste(object_loc, 'combined.v3.20210629.ct.rds', sep = '')

# load objects
v2_combined <- readRDS(v2_combined_loc)
v3_combined <- readRDS(v3_combined_loc)

# add pathways
v2_combined <- add_pathway_module_score(v2_combined, pathway_annotation_il6, assay='SCT', name='IL6')
v2_combined <- add_pathway_module_score(v2_combined, pathway_annotation_il6, assay='SCT', name='IL6')

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

