#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_simplify_limix_covariate_files.R
# Function: simplify the covariates files that are used as input for LIMIX
############################################################################################################################

####################
# libraries        #
####################

####################
# Functions        #
####################

#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param covariate_files_folder where the current covariate files are
#' @param columns_to_keep columns in the covariate files to keep after filtering
#' @param id_column the column to use as the ID column that is matched to expression files
#' @param regex_pattern pattern for finding the files in the output directory
#' @param name_pattern for the name of the new file, what to replace
#' @param name_replacement for the name of the new file, what to replace with
#' @returns 0 if succesfull
#' filter_covariates_files(qtl_input_low_loc, columns_to_keep, id_column)
filter_covariates_files <- function(covariate_files_folder, columns_to_keep, id_column, regex_pattern='*covariates*', name_pattern='covariates', name_replacement='covariates_simplified', make_categorical_values_binary=T, categorial_exceptions=c('Donor_Pool'), unknown_name='unknown') {
  # get the current covariate files
  covariate_file_locs <- list.files(covariate_files_folder, pattern = regex_pattern, full.names = F)
  # check each file
  for (covariate_file_loc in covariate_file_locs) {
    # read the original file
    covariates <- read.table(paste(covariate_files_folder, covariate_file_loc, sep = ''), header = T, check.names = F, comment.char = '', stringsAsFactors = F)
    # keep only the ones we care about
    covariates <- covariates[, c(id_column, setdiff(columns_to_keep, id_column))]
    # remove duplicates
    covariates <- unique(covariates)
    # do categorical to binary mapping if requested
    if (make_categorical_values_binary) {
      # save which columns where categorical
      categoricals <- c()
      # get the original columns
      columns_original <- colnames(covariates)
      # check each column
      for (column in columns_original) {
        # only do this if it is not in the exception list
        if (!(column %in% categorial_exceptions)) {
          # get the type of the column
          if (typeof(covariates[[column]]) == 'character') {
            # turn NA into 'unknown'
            covariates[is.na(covariates[[column]]), column] <- unknown_name
            # add to the categorical list
            categoricals <- c(categoricals, column)
            # get the unique combinations in the column
            categories <- unique(covariates[[column]])
            # if there is only one category, we cannot do anything
            if (length(categories) == 1) {
              message(paste('only one category for', column, 'will be dropped'))
            }
            else{
              # we need all the categories minus one, as it if is not one of the others, it automatically is the one we didn't include
              categories <- categories[1:(length(categories) - 1)]
              # now check for each category
              for (category in categories) {
                covariates[[paste(column, category, sep = '_')]] <- ifelse(covariates[[column]] == category, 1, 0)
              }
            }
          }
        }
      }
      # remove original columns
      covariates[categoricals] <- NULL
    }
    # create a new output file
    filtered_covariate_file_loc <- gsub(name_pattern, name_replacement, covariate_file_loc)
    # the new output file loc
    filtered_covariate_file_path <- paste(covariate_files_folder, filtered_covariate_file_loc, sep = '')
    # write the result
    write.table(covariates, filtered_covariate_file_path, sep = '\t', row.names = F, col.names = T, quote = F)
  }
  return(0)
}

####################
# Main Code        #
####################

# the columns we want to keep
columns_to_keep <- c('Donor_Pool', 'sample_condition', paste('PC', 1:10, sep = ''))
# the column to have as the first column, that will be the sample designation
id_column <- 'Donor_Pool'

# location of the inputs
qtl_input_loc <- '/scratch/p287578/projects/venema-2022/ongoing/qtl/eqtl/sc-eqtlgen/input/elmentaite_adult_martin_immune/'
qtl_input_low_loc <- paste(qtl_input_loc, 'cell_type_low_inflammationsplit/', sep = '')
qtl_input_med_loc <- paste(qtl_input_loc, 'cell_type_med_inflammationsplit/', sep = '')
qtl_input_medhigh_loc <- paste(qtl_input_loc, 'cell_type_medhigh_inflammationsplit/', sep = '')
qtl_input_high_loc <- paste(qtl_input_loc, 'cell_type_safe_inflammationsplit/', sep = '')

# filter covariates
filter_covariates_files(qtl_input_low_loc, columns_to_keep, id_column, regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_med_loc, columns_to_keep, id_column, regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_medhigh_loc, columns_to_keep, id_column, regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_high_loc, columns_to_keep, id_column, regex_pattern='.covariates.txt')

# the columns we want to keep
columns_to_keep <- c('Donor_Pool', 'sample_condition')
# filter covariates
filter_covariates_files(qtl_input_low_loc, columns_to_keep, id_column, name_replacement='covariates_status', regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_med_loc, columns_to_keep, id_column, name_replacement='covariates_status', regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_medhigh_loc, columns_to_keep, id_column, name_replacement='covariates_status', regex_pattern='.covariates.txt')
filter_covariates_files(qtl_input_high_loc, columns_to_keep, id_column, name_replacement='covariates_status', regex_pattern='.covariates.txt')
