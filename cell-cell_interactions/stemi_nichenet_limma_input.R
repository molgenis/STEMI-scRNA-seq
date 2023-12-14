#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name:  stemi_nichenet_limma_input.R
# Function: create limma output files for each celltype and condition combination, to use for nichenet
############################################################################################################################

####################
# libraries        #
####################

# required for nichenet
library(Seurat)
# limma DE dependencies
library(variancePartition)
library(edgeR)
library(BiocParallel)
# convert count matrices
library(Matrix)


####################
# Functions        #
####################


#' create a column which combines the values of multiple other columns
#' 
#' @param dataframe dataframe that has the original column, to which to add the aggregate column
#' @param columns the columns to combine into a new column
#' @returns the original dataframe, with a new column 'all_aggregates', which is the combination of the supplied columns
#' new_df <- combine_columns(df, c('V1', 'V3'))
combine_columns <- function(dataframe, columns) {
  # check each column
  for (column in columns) {
    # check if we were already aggregating
    if ('all_aggregates' %in% colnames(dataframe)) {
      # add the new column
      dataframe[['all_aggregates']] <- paste(dataframe[['all_aggregates']], dataframe[[column]], sep = '_')
    }
    # otherwise we need to start our aggregation
    else{
      dataframe[['all_aggregates']] <- dataframe[[column]]
    }
  }
  return(dataframe)
}


#' get the number of cells describing a combination of values from the seurat metadata
#' 
#' @param aggregate_df the dataframe with the combinations that were aggregated over
#' @param aggregated_columns, the columns that were aggregated over
#' @param seurat_metadata the metadata of a seurat object from which to get cell numbers for the combinations
#' @returns the original dataframe, with a 'nr' column, denoting the number of cells for that combination
#' cell_numbers <- get_nr_cells_aggregate_combination(aggretates, pbmc_meta.data)
get_nr_cells_aggregate_combination <- function(aggregate_df, aggregate_columns, seurat_metadata) {
  # now use table to get the number of entries of these aggregates in the metadata
  aggregate_numbers <- data.frame(table(seurat_metadata[, aggregate_columns]))
  # add a new column for the combination of the aggregates
  aggregate_df <- combine_columns(aggregate_df, aggregate_columns)
  # in our numbers as well
  aggregate_numbers <- combine_columns(aggregate_numbers, aggregate_columns)
  # now join the frequencies from the aggragate numbers onto the aggregate columns
  aggregate_df[['nr']] <- aggregate_numbers[match(aggregate_df[['all_aggregates']], aggregate_numbers[['all_aggregates']]), 'Freq']
  # now we are sure that the cell numbers are in the same order as the aggregated metadata
  # remove the column we created
  aggregate_df[['all_aggregates']] <- NULL
  # and return the result
  return(aggregate_df)
}

get_aggregate_object <- function(seurat_object, aggregate_columns, assay_oi='SCT') {
  # add aggregation columns
  for(aggregate_column in aggregate_columns) {
    if ('COMBINATION' %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[['COMBINATION']] <- paste(seurat_object@meta.data[['COMBINATION']], seurat_object@meta.data[[aggregate_column]], sep = '_SEPARATOR_')
    }
    else {
      seurat_object@meta.data[['COMBINATION']] <- seurat_object@meta.data[[aggregate_column]]
    }
  }
  # go the aggregation
  seurat_agg <- AggregateExpression(seurat_object, assays = assay_oi, return.seurat = T, group.by = c('COMBINATION'), slot = 'counts')
  # split back the combined metadata columns
  seurat_agg@meta.data[aggregate_columns] <- str_split_fixed(rownames(seurat_agg@meta.data), '_SEPARATOR_', length(aggregate_columns))
  # and clean up the names
  seurat_agg <- RenameCells(seurat_agg, new.names = gsub('_SEPARATOR_','_', colnames(seurat_agg)))
  # next get the cell numbers for each observation, these will follow the order of the original aggregated metadata
  cell_numbers <- get_nr_cells_aggregate_combination(seurat_agg@meta.data[, aggregate_columns], aggregate_columns, seurat_object@meta.data)
  # add to the metadata
  seurat_agg[['ncells']] <- cell_numbers[['nr']]
  return(seurat_agg)
}


do_dream <- function(aggregated_seurat_object, contrast, condition_1, condition_2, output_loc, fixed_effects=c(), random_effects=c(), minimal_cells=10, minimal_complexity=10000, verbose=T) {
  # filter by the number of cells if requested
  if (minimal_cells > 0) {
    # get the indices of where the cell numbers are above this
    indices_above_threshold <- which(aggregated_seurat_object@meta.data[['ncells']] >= minimal_cells)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(aggregated_seurat_object@meta.data)), 'entries, ', length(indices_above_threshold), 'contained more cells than the', as.character(minimal_cells), 'threshold'))
    }
    # do the actual filtering
    aggregated_seurat_object <- aggregated_seurat_object[, indices_above_threshold]
  }
  # filter by the complexity if requested
  if (minimal_complexity > 0) {
    # calculate the complexity first
    complexity <- data.frame(colSums(aggregated_seurat_object@assays[[assay_oi]]@counts))
    # set better column names
    colnames(complexity) <- c('complexity')
    # get the indices of where the complexity is above the threshold
    indices_above_complexity <- which(complexity[['complexity']] >= minimal_complexity)
    # report how many
    if (verbose) {
      message(paste('of', as.character(nrow(complexity)), 'entries, ', length(indices_above_complexity), 'contained more UMIs than the', as.character(minimal_complexity), 'threshold'))
    }
    # and filter
    aggregated_seurat_object <- aggregated_seurat_object[, indices_above_complexity]
  }
  # filter genes by number of counts
  isexpr = rowSums(cpm(aggregated_seurat_object@assays$SCT@counts)>0.1) >= 5
  
  # Standard usage of limma/voom
  geneExpr = DGEList( aggregated_seurat_object@assays$SCT@counts[isexpr,] )
  geneExpr = calcNormFactors( geneExpr )
  
  # Specify parallel processing parameters
  # this is used implicitly by dream() to run in parallel
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  # if there are no fixed effects, we'll just use the contrast
  if (length(fixed_effects) == 0) {
    fixed_effects <- c(contrast)
  }
  
  # paste together the model
  model_formula <- '~ 0 '
  for(fixed_effect in fixed_effects){
    model_formula <- paste(model_formula, fixed_effect, sep = ' + ')
  }
  for(random_effect in random_effects){
    model_formula <- paste(model_formula, ' + (1|', random_effect, ')', sep = '')
  }
  
  if(verbose){
    print(paste('formula:', model_formula))
  }
  # and turn into a formula
  form <- as.formula(model_formula)
  tryCatch(
    {
      # estimate weights using linear mixed model of dream
      vobjDream = voomWithDreamWeights( counts = geneExpr, formula = form, data = aggregated_seurat_object@meta.data, weights = aggregated_seurat_object@meta.data[['ncells']] ) # the cell numbers are in the same order as the metadata, and as such can be passed like this
      tryCatch(
        {
          # define and then cbind contrasts
          L = getContrast( vobjDream, form, aggregated_seurat_object@meta.data, c(paste(contrast, condition_1, sep=''), paste(contrast, condition_2, sep='')))
          # fit contrast
          fit = dream( vobjDream, form, aggregated_seurat_object@meta.data, L)
          
          # if there are no random effects, we need to use ebayes
          if (is.null(random_effects) | length(random_effects) == 0) {
            message('no random effects supplied, running ebayes')
            fit <- eBayes(fit)
          }
          
          # grab the exact fit
          limma_result <- topTable(fit, coef='L1', number=length(fit$F.p.value))
          
          # add bonferroni adjustment
          limma_result[['p.bonferroni']] <- p.adjust(limma_result[['P.Value']])
          
          # add some statistics
          result_stats_list <- list()
          # check each condition
          result_stats_list[['combination']] <- data.frame(combination=rep(paste(contrast, paste(c(condition_1, condition_2), collapse='-'), sep = '.'), times = nrow(limma_result)))
          for (condition in c(condition_1, condition_2)) {
            # get the cell numbers for the condition
            cell_numbers_condition <- aggregated_seurat_object@meta.data[aggregated_seurat_object@meta.data[[contrast]] == condition, ]
            result_stats_list[[paste('ndonor', condition, sep = '_')]] <- data.frame(ndonor=rep(nrow(cell_numbers_condition), times = nrow(limma_result)))
            result_stats_list[[paste('ncell', condition, sep = '_')]] <- data.frame(ncells=rep(
              paste(as.character(min(cell_numbers_condition[['ncells']])),
                    as.character(quantile(cell_numbers_condition[['ncells']])[['25%']]),
                    as.character(quantile(cell_numbers_condition[['ncells']])[['50%']]),
                    as.character(quantile(cell_numbers_condition[['ncells']])[['75%']]),
                    as.character(max(cell_numbers_condition[['ncells']])),
                    sep = ';'
              ), times = nrow(limma_result)))
          }
          # merge all
          result_stats <- do.call('cbind', result_stats_list)
          colnames(result_stats) <- names(result_stats_list)
          
          # add combination as first column
          limma_result <- cbind(result_stats, limma_result)
          
          # set an output location
          limma_output_loc <- paste(output_loc, contrast, '.', paste(c(condition_1, condition_2), collapse = '_'), '.tsv.gz', sep = '')
          
          # also write the model we used
          limma_formula_loc <- paste(output_loc, contrast, '.', paste(c(condition_1, condition_2), collapse = '_'), '.formula', sep = '')
          
          if(verbose){
            print(paste('writing result', limma_output_loc))
          }
          
          # write the result
          write.table(limma_result, gzfile(limma_output_loc), sep = '\t', row.names = T)
          # and the formula
          write.table(model_formula, limma_formula_loc, row.names = F, col.names = F)
        }, error=function(cond) {
          print(paste('analysis failed in', paste(c(condition_1, condition_2))))
          message(cond)
        }
      )
    }
  )
}


####################
# Main Code        #
####################

# location of the object
objects_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/'
cardio_integrated_loc <- paste(objects_loc, 'cardio.integrated.20210301.rds', sep = '')

# read the object
cardio_integrated <- readRDS(cardio_integrated_loc)
cardio_integrated <- UpdateSeuratObject(cardio_integrated)
DefaultAssay(cardio_integrated) <- 'SCT'

# convert the object into a pseudobulk
cardio_agg <- get_aggregate_object(cardio_integrated, c('chem', 'timepoint.final', 'assignment.final', 'cell_type_lowerres', 'age', 'gender'))
# make age a number again
cardio_agg@meta.data[['age']] <- as.numeric(cardio_agg@meta.data[['age']])
# add a column that is the combination of the celltype and the condition
cardio_agg@meta.data[['celltype_condition']] <- paste(cardio_agg@meta.data[['cell_type_lowerres']], cardio_agg@meta.data[['timepoint.final']], sep = '_')
# keep only the major cell types
cardio_agg <- cardio_agg[, cardio_agg@meta.data[['cell_type_lowerres']] %in% c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')]
# do each combination
all_combs <- unique(cardio_agg@meta.data[['celltype_condition']])
comb_length <- length(all_combs)
for (i in 1:(comb_length - 1)) {
  for (i2 in (i+1):comb_length) {
    print(paste(c(all_combs[i], all_combs[i2])))
    do_dream(cardio_agg[, cardio_agg@meta.data$celltype_condition %in% c(all_combs[i], all_combs[i2]), ], 'celltype_condition', all_combs[i2], all_combs[i], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_cell_interactions/nichenet/limma_de/', fixed_effects=c('celltype_condition', 'chem', 'age', 'gender'), random_effects = c())
  }
}
