######################
# libraries          #
######################

library(Seurat)
library(Matrix)

####################
# Functions        #
####################

create_features_files_bulk <- function(seurat_object,  output_loc, symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", assay = "RNA", prepend_1 = T){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # get the individuals in the object
  individuals <- unique(seurat_object@meta.data[,sample.id.column.name])
  # remove any na individuals (sometimes caused if using demux assignments)
  individuals <- individuals[!is.na(individuals)]
  
  # create a matrix where we will store the mean expressions
  mean_expression_matrix <- matrix(nrow=nrow(seurat_object), ncol = length(individuals), 
                                   dimnames = list(rownames(seurat_object), individuals))
  # go through the individuals
  for (individual in individuals) {
    # take into account that there might be zero expression
    if (sum(is.na(seurat_object@meta.data[,sample.id.column.name]) | seurat_object@meta.data[,sample.id.column.name] == individual) == 0) {
      mean_expression_matrix[,individual] <- 0
    }
    else {
      # grab the cells of this cell type and this individual
      cells_individual <- seurat_object[,!(is.na(seurat_object@meta.data[,sample.id.column.name])) & seurat_object@meta.data[,sample.id.column.name] == individual]
      # calculate the mean expression for each gene
      mean_expression_matrix[,individual] <- rowMeans(cells_individual)
    }
  }
  # remap symbols to ensg numbers if necessary
  if (symbols.to.ensg) {
    genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
    genes$V2 <- gsub("_", "-", make.unique(genes$V2))
    rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
  }
  if(prepend_1){
    # prepend the '1_' to the id 
    colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
  }
  
  # write our table of means for this cell type
  write.table(mean_expression_matrix, 
              file = paste0(output_loc, "bulk", "_expression", ".tsv"),
              quote = F, sep = "\t", col.names = NA)
}

create_feature_files_per_condition_bulk <- function(seurat_object, output_loc, condition.column.name = "timepoint", symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", assay = "RNA", cell_type_column = "cell_type", prepend_1 = T){
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    # do the feature file creation for this condition
    create_features_files_bulk(cells_condition, output_loc_condition, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping = symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, assay=assay, prepend_1 = prepend_1)
  }
}

create_features_files <- function(seurat_object,  output_loc, cell_types_to_output = NULL, symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment", cell_type_column = "cell_type", assay = "RNA", prepend_1 = T){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # get the individuals in the object
  individuals <- unique(seurat_object@meta.data[,sample.id.column.name])
  # remove any na individuals (sometimes caused if using demux assignments)
  individuals <- individuals[!is.na(individuals)]
  # if the user did not specify the cell types to output, just do all of them
  if(is.null(cell_types_to_output)){
    cell_types_to_output = unique(seurat_object@meta.data[[cell_type_column]])
  }
  # we will want to report on the number of cells per cell type, so that needs to be stored in a matrix
  cell_counts <- matrix(nrow=length(cell_types_to_output), ncol = length(individuals), 
                        dimnames = list(cell_types_to_output, individuals))
  # go through the cell types we want to output
  for (celltype in cell_types_to_output) {
    # grab the cells specific to that cell type
    cells_cell_type <- seurat_object[,seurat_object@meta.data[cell_type_column] == celltype]
    # create a matrix where we will store the mean expressions
    mean_expression_matrix <- matrix(nrow=nrow(cells_cell_type), ncol = length(individuals), 
                                     dimnames = list(rownames(cells_cell_type), individuals))
    # go through the individuals
    for (individual in individuals) {
      # take into account that there might be zero expression
      if (sum(is.na(cells_cell_type@meta.data[,sample.id.column.name]) | cells_cell_type@meta.data[,sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[,individual] <- 0
        cell_counts[celltype,individual] <- 0
      }
      else {
        # grab the cells of this cell type and this individual
        cells_cell_type_individual <- cells_cell_type[,cells_cell_type@meta.data[,sample.id.column.name] == individual]
        # calculate the mean expression for each gene
        mean_expression_matrix[,individual] <- rowMeans(cells_cell_type_individual)
        # get the number of cells of this cell type for this individual
        cell_counts[celltype,individual] <- ncol(cells_cell_type_individual)
      }
    }
    # remap symbols to ensg numbers if necessary
    if (symbols.to.ensg) {
      genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
      genes$V2 <- gsub("_", "-", make.unique(genes$V2))
      rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
    }
    if(prepend_1){
      # prepend the '1_' to the id 
      colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
    }
    # write our table of means for this cell type
    write.table(mean_expression_matrix, 
                file = paste0(output_loc, celltype, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)
  }
  # write out table of 
  write.table(cell_counts, 
              file = paste0(output_loc, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

create_feature_files_per_condition <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint", sample.id.column.name="assignment", cell_type_column = "cell_type", assay = "RNA", symbols.to.ensg = F, symbols.to.ensg.mapping="genes.tsv", prepend_1 = T){
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # unless we wnat only specific conditions, we'll just to those
  if(!(is.null(conditions_to_output))){
    conditions <- conditions_to_output
  }
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    print(paste("grabbing condition", condition))
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    print(paste("finished grabbing condition", condition))
    # do the feature file creation for this condition
    create_features_files(cells_condition, cell_types_to_output = NULL, output_loc = output_loc_condition, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping=symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = prepend_1)
  }
}

######################
# main code          #
######################

#locations of objects
objects_loc <- "/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
cardio.integrated_loc <- paste(objects_loc, "cardio_integrated_20pcs_ctd_final.rds", sep = '')

# locations of the features
#features_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/'
features_loc <- '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/'
features_loc_stemi_v2_full_res <- paste(features_loc, 'stemi_v2_full_res/', sep = '')
features_loc_stemi_v3_full_res <- paste(features_loc, 'stemi_v3_full_res/', sep = '')
features_loc_stemi_v2_lowerres <- paste(features_loc, 'stemi_v2_lowerres/', sep = '')
features_loc_stemi_v3_lowerres <- paste(features_loc, 'stemi_v3_lowerres/', sep = '')
features_loc_stemi_v2_lowestres <- paste(features_loc, 'stemi_v2_lowestres/', sep = '')
features_loc_stemi_v3_lowestres <- paste(features_loc, 'stemi_v3_lowestres/', sep = '')

# for the inhouse eQTL-mapping pipeline, we currently have ENSG numbers instead of gene numbers
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"

# read the object
cardio.integrated <- readRDS(cardio.integrated_loc)

# subset stemi v2
stemi_v2 <- subset(cardio.integrated, subset = orig.ident == 'stemi_v2')
DefaultAssay(stemi_v2) <- 'SCT'
# create the features files
create_feature_files_per_condition_bulk(stemi_v2, features_loc_stemi_v2_full_res, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type", prepend_1 = F)
create_feature_files_per_condition(stemi_v2, features_loc_stemi_v2_full_res, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)
create_feature_files_per_condition_bulk(stemi_v2, features_loc_stemi_v2_lowerres, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type_lowerres", prepend_1 = F)
create_feature_files_per_condition(stemi_v2, features_loc_stemi_v2_lowerres, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)
create_feature_files_per_condition_bulk(stemi_v2, features_loc_stemi_v2_lowestres, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type_lowestres", prepend_1 = F)
create_feature_files_per_condition(stemi_v2, features_loc_stemi_v2_lowestres, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type_lowestres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)

# subset stemi v3
stemi_v3 <- subset(cardio.integrated, subset = orig.ident == 'stemi_v3')
DefaultAssay(stemi_v3) <- 'SCT'
# create the features files
create_feature_files_per_condition_bulk(stemi_v3, features_loc_stemi_v3_full_res, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type", prepend_1 = F)
create_feature_files_per_condition(stemi_v3, features_loc_stemi_v3_full_res, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)
create_feature_files_per_condition_bulk(stemi_v3, features_loc_stemi_v3_lowerres, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type_lowerres", prepend_1 = F)
create_feature_files_per_condition(stemi_v3, features_loc_stemi_v3_lowerres, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type_lowerres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)
create_feature_files_per_condition_bulk(stemi_v3, features_loc_stemi_v3_lowestres, condition.column.name = "timepoint.final", symbols.to.ensg=T, symbols.to.ensg.mapping = gene_to_ens_mapping, sample.id.column.name="assignment.final", assay = "SCT", cell_type_column = "cell_type_lowestres", prepend_1 = F)
create_feature_files_per_condition(stemi_v3, features_loc_stemi_v3_lowestres, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type_lowestres", assay = "SCT", symbols.to.ensg = T, symbols.to.ensg.mapping=gene_to_ens_mapping, prepend_1 = F)

