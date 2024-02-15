#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_create_celltype_composition_tables.R
# Function: create tables used for cell type composition analysis
############################################################################################################################

####################
# libraries        #
####################

library(Seurat)

####################
# Functions        #
####################

metadata_to_composition_table <- function(seurat_metadata, cell_type_column='cell_type_final', split_columns=c('donor_final', 'inflammation_status')){
  # replace the cell types with something that is safe
  seurat_metadata[[cell_type_column]] <- gsub(' |/', '_', seurat_metadata[[cell_type_column]])
  seurat_metadata[[cell_type_column]] <- gsub('-', '_negative', seurat_metadata[[cell_type_column]])
  seurat_metadata[[cell_type_column]] <- gsub('\\+', '_positive', seurat_metadata[[cell_type_column]])
  seurat_metadata[[cell_type_column]] <- gsub('\\)', '', seurat_metadata[[cell_type_column]])
  seurat_metadata[[cell_type_column]] <- gsub('\\(', '', seurat_metadata[[cell_type_column]])
  # if we have no split columns, we can get the table
  if(length(split_columns) == 0){
    # use tabel to get the numbers
    numbers <- t(table(seurat_metadata[, c(cell_type_column), drop = F]))
    # extract the cell names
    cell_names <- colnames(numbers)
    # extract the raw values
    values_vector <- as.vector(numbers)
    # make a single row dataframe
    numbers_table <- data.frame(matrix(, nrow = 1, ncol = length(values_vector), dimnames = list(NA, cell_names)))
    # set the values
    numbers_table[1, cell_names] <- values_vector
    return(numbers_table)
  }
  else{
    # if there is a split column left, grab the last one
    column_to_split_on <- split_columns[length(split_columns)]
    # we will save a result for this iteration
    numbers_table <- NULL
    # then check each variant
    for(column_category in unique(seurat_metadata[[column_to_split_on]])){
      # subset to that category
      seurat_metadata_category <- seurat_metadata[seurat_metadata[[column_to_split_on]] == column_category, ]
      # recursively call this function, with one less split column (as we just subsetted for that one)
      counts_category <- metadata_to_composition_table(seurat_metadata_category, cell_type_column, setdiff(split_columns, c(column_to_split_on)))
      # add that these counts are for what we just subsetted for
      counts_category[[column_to_split_on]] <- column_category
      # now add result with the rest of the categories of this level
      if(is.null(numbers_table)){
        numbers_table <- counts_category
      }
      else{
        # we have to make sure that we have all cell types in all categories
        missing_in_category <- setdiff(colnames(numbers_table), colnames(counts_category))
        missing_in_other_categories <- setdiff(colnames(counts_category), colnames(numbers_table))
        # we will add these as '0' in both
        counts_category[, missing_in_category] <- 0
        numbers_table[, missing_in_other_categories] <- 0
        # now they can be safely rbound
        numbers_table <- rbind(numbers_table, counts_category)
      }
    }
    return(numbers_table)
  }
}


add_lower_classification <- function(seurat_object, reclassification_mapping, mapping_original_column, mapping_reclass_column, metadata_original_column, metadata_reclassification_column){
  # add the new column
  seurat_object@meta.data[[metadata_reclassification_column]] <- NA
  # get each cell type in the data
  metadata_original_cts <- unique(seurat_object@meta.data[[metadata_original_column]])
  # and the originals in the mapping
  reclassification_original_cts <- unique(reclassification_mapping[[mapping_original_column]])
  # we can only map what is present in both
  originals_both <- intersect(metadata_original_cts, reclassification_original_cts)
  # check what is missing
  only_metadata <- setdiff(metadata_original_cts, reclassification_original_cts)
  only_mapping <- setdiff(reclassification_original_cts, metadata_original_cts)
  # warn what is missing
  if(length(only_metadata) > 0){
    print('some celltypes only in metadata')
    print(only_metadata)
  }
  if(length(only_mapping) > 0){
    print('some celltypes only in remapping ')
    print(only_mapping)
  }
  # check each cell type
  for(celltype_original in originals_both){
    # get the appropriate remapping
    celltype_remapped <- reclassification_mapping[reclassification_mapping[[mapping_original_column]] == celltype_original, mapping_reclass_column]
    # now remap in the metadata
    seurat_object@meta.data[seurat_object@meta.data[[metadata_original_column]] == celltype_original, metadata_reclassification_column] <- celltype_remapped
  }
  return(seurat_object)
}


make_celltypes_safe <- function(cell_types){
  # get a safe file name
  cell_type_safes <- gsub(' |/', '_', cell_types)
  cell_type_safes <- gsub('-', '_negative', cell_type_safes)
  cell_type_safes <- gsub('\\+', '_positive', cell_type_safes)
  cell_type_safes <- gsub('\\)', '', cell_type_safes)
  cell_type_safes <- gsub('\\(', '', cell_type_safes)
  return(cell_type_safes)
}


####################
# Main Code        #
####################


# location of objects
objects_loc <- '/groups/umcg-franke-scrna/tmp03/releases/blokland-2020/v1/seurat/'
stemi_object_loc <- paste(objects_loc, 'cardio.integrated.20210301.rds', sep = '')
# read the Seurat object
stemi <- readRDS(stemi_object_loc)

# location of the cell number tables
cell_number_tables_loc <- '/groups/umcg-franke-scrna/tmp03/releases/blokland-2020/v1/cell_type_composition/cell_number_tables/'
cell_number_tables_highres_loc <- paste(cell_number_tables_loc, 'stemi_cell_numbers_highres.tsv', sep = '')
cell_number_tables_lowres_loc <- paste(cell_number_tables_loc, 'stemi_cell_numbers_lowres.tsv', sep = '')

# for the final cell type classification
stemi_cell_numbers_ <- metadata_to_composition_table(stemi@meta.data, cell_type_column = 'cell_type', split_columns = c('timepoint.final', 'assignment.final'))
# add the sex and age as well from the cell data. Since it's the same for each cell of the same donor, we can just use the first match
stemi_cell_numbers_[, c('age', 'sex')] <- stemi@meta.data[match(stemi_cell_numbers_[['assignment.final']], stemi@meta.data[['assignment.final']]), c('age', 'gender')]
# remove empty entries
stemi_cell_numbers_ <- stemi_cell_numbers_[apply(stemi_cell_numbers_, 1, function(x){sum(is.na(x))}) == 0, ]
write.table(stemi_cell_numbers_, cell_number_tables_highres_loc, sep = '\t', row.names = F, col.names = T, quote = F)

# for the final cell type classification
stemi_cell_numbers_lowres <- metadata_to_composition_table(stemi@meta.data, cell_type_column = 'cell_type_lowerres', split_columns = c('timepoint.final', 'assignment.final'))
# add the sex and age as well from the cell data. Since it's the same for each cell of the same donor, we can just use the first match
stemi_cell_numbers_lowres[, c('age', 'sex')] <- stemi@meta.data[match(stemi_cell_numbers_lowres[['assignment.final']], stemi@meta.data[['assignment.final']]), c('age', 'gender')]
# remove empty entries
stemi_cell_numbers_lowres <- stemi_cell_numbers_lowres[apply(stemi_cell_numbers_lowres, 1, function(x){sum(is.na(x))}) == 0, ]
write.table(stemi_cell_numbers_lowres, cell_number_tables_lowres_loc, sep = '\t', row.names = F, col.names = T, quote = F)

