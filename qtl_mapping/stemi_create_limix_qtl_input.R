#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Marc-Jan Bonder
# Name: stemi_create_limix_qtl_input.R
# Function: create the limix-QTL compatible input files from the Seurat object
############################################################################################################################

####################
# libraries        #
####################

library(data.table)
library(Seurat)
library(matrixStats)
library(textTinyR) # NOT IN CONTAINER
library(pbapply)

####################
# Functions        #
####################


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object_metadata the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param pool_column the seurat metadata column that denotes the pool that the samples were run in
#' @param condition_column the seurat metadata column that denotes the condition of the sample (inflamed, non-inflamed)
#' @param celltype_column rthe seurat metadata column that denotes the celltype of the cell
#' @returns a list per cell type, with for each cell type a metadata summary table
#' metadata_per_celltype <- create_metadata_tables(seurat_object@meta.data, donor_annotation_psam, participant_column = 'soup_final_sample_assignment')
create_metadata_tables <- function(seurat_object_metadata, psam, participant_column='donor_final', pool_column='day', condition_column='inflammation_status', celltype_column='cell_type_safe', join_pools=T) {
  # we need a separate psam for each cell type
  psam_per_celltype <- list()
  # add the sequencing platform to the psam
  psam[['sequencing_platform']] <- 'BGIseq'
  # some variables are more complicated, as they are a combination of multiple variables
  psam[['sequencing_run']] <- NA # lane
  psam[['sequencing_lane']] <- NA # lane
  psam[['scrna_platform']] <- '10x_v2.0'
  psam[['plate_base']] <- 'N'
  psam[['umi_based']] <- 'Y'
  psam[['biomaterial']] <- 'PBMC'
  psam[['sorting']] <- 'N'
  psam[['cell_treatment']] <- 'UT'
  psam[['sample_condition']] <- NA # inflamed/non-inflamed
  psam[['Donor_Pool']] <- NA # participant;;lane
  # now check the values that vary between individuals
  for (row_i in 1:nrow(psam)) {
    # extract the participant
    participant <- psam[row_i, 'IID']
    # check which lanes this participant is in
    lanes <- unique(seurat_object_metadata[seurat_object_metadata[[participant_column]] == participant, pool_column])
    # check which conditions this sample was in
    conditions <- unique(seurat_object_metadata[seurat_object_metadata[[participant_column]] == participant, condition_column])
    # check if the lanes are not empty
    if (length(lanes) == 0) {
      # set as unknown if this is the case
      lanes <- c(NA)
    }
    # and for conditions
    if (length(conditions) == 0) {
      conditions <- c(NA)
    }
    # now set these values
    psam[row_i, 'sequencing_run'] <- paste(lanes, collapse=',')
    psam[row_i, 'sequencing_lane'] <- paste(lanes, collapse = ',')
    psam[row_i, 'sample_condition'] <- paste(conditions, collapse = ',')
    psam[row_i, 'Donor_Pool'] <- paste(participant, lanes, sep=';;', collapse = ',')
  }
  if (!join_pools) {
    # if we shouldn't join pools, let's get those specific ones
    indices_multiple_pools <- grep(',', psam[['Donor_Pool']])
    rows_multiple_pools <- psam[indices_multiple_pools, ]
    # and remove them from the original
    rows_single_pools <- psam[-indices_multiple_pools, ]
    # create a new list
    rows_split <- list()
    for (row_i in 1:nrow(rows_multiple_pools)){
      # extract the donor pools
      donor_pools <- strsplit(rows_multiple_pools[row_i, 'Donor_Pool'], ',')[[1]]
      # and the sample condition
      sample_conditions <- strsplit(rows_multiple_pools[row_i, 'sample_condition'], ',')[[1]]
      # and the sequencing lane
      sequencing_lanes <- strsplit(rows_multiple_pools[row_i, 'sequencing_lane'], ',')[[1]]
      # and the sequencing run
      sequencing_runs <- strsplit(rows_multiple_pools[row_i, 'sequencing_run'], ',')[[1]]
      # now add those as split column
      for (i in 1:length(donor_pools)) {
        # get individual pool and sample condition
        donor_pool <- donor_pools[i]
        condition <- sample_conditions[i]
        lane <- sequencing_lanes[i]
        run <- sequencing_runs[i]
        # create a new row
        new_row <- rows_multiple_pools[row_i, ]
        new_row[, 'Donor_Pool'] <- donor_pool
        new_row[, 'sample_condition'] <- condition
        new_row[, 'sequencing_lane'] <- lane
        new_row[, 'sequencing_run'] <- run
        # add it to the new split rows
        rows_split[[donor_pool]] <- new_row
      }
    }
    # merge the new rows
    df_rows_split <- do.call(rbind, rows_split)
    # then merge the old and the new
    psam <- rbind(rows_single_pools, df_rows_split)
  }
  # now check each cell type
  for (cell_type in unique(seurat_object_metadata[[celltype_column]])) {
    # subset the metadata to that celltype
    metadata_celltype <- seurat_object_metadata[seurat_object_metadata[[celltype_column]] == cell_type, ]
    # create a new psam for this cell type
    psam_celltype <- psam
    # if we joined the pools, we can just use the donor to get the number of cells
    if (join_pools) {
      # now use table to get the cell type numbers
      cell_numbers_per_donor <- data.frame(table(metadata_celltype[[participant_column]]))
      # now add the cell counts
      psam_celltype[['CellCount']] <- cell_numbers_per_donor[match(psam_celltype[['IID']], cell_numbers_per_donor[['Var1']]), 'Freq']
    }
    # if we have the donors entered multiple times with different pools, we need to get the cells for each specific combination
    else{
      # paste together the donor and the pool
      metadata_celltype[['Donor_Pool']] <- paste(metadata_celltype[[participant_column]], metadata_celltype[[pool_column]], sep = ';;')
      # and instead get the cell numbers like that
      cell_numbers_per_donor_pool <- data.frame(table(metadata_celltype[['Donor_Pool']]))
      # and add that
      psam_celltype[['CellCount']] <- cell_numbers_per_donor_pool[match(psam_celltype[['Donor_Pool']], cell_numbers_per_donor_pool[['Var1']]), 'Freq']
    }
    # turn NA into zero
    psam_celltype[is.na(psam_celltype[['CellCount']]), 'CellCount'] <- 0
    # add to list of psams
    psam_per_celltype[[cell_type]] <- psam_celltype
  }
  return(psam_per_celltype)
}


#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object the metadata (from Seurat) to use to create a metadata annotation file
#' @param participant_column the seurat metadata column that denotes the participant
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param batch_column the batch the sample was processed in (optional)
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param verbose print progress or not
#' @returns a list per cell type, each cell type has a list with the raw pseudobulk expression, the filtered pseudobulk expression, and the pcs
#' expression_per_celltype <- create_aggregated_expression_matrices(seurat_object, participant_column = 'soup_final_sample_assignment')
create_aggregated_expression_matrices <- function(seurat_object, participant_column='donor_final', celltype_column='cell_type_safe', batch_column=NULL, min_cell_number=5, min_numi=200, npcs=10, verbose=T) {
  # subset object if min_numi parameter is given
  if (!is.null(min_numi) & !is.na(min_numi) & min_numi > 0) {
    seurat_object <- seurat_object[, seurat_object@meta.data[['nFeature_SCT']] >= min_numi]
  }
  # if the batch column is non-empty, we need to make a new column that is the join of the batch column and the participant ID
  if (!is.null(batch_column)) {
    seurat_object@meta.data[['participant_pool']] <- paste(seurat_object@meta.data[[participant_column]], seurat_object@meta.data[[batch_column]], sep = ';;')
    participant_column <- 'participant_pool'
  }
  # we'll store the aggregated count matrices per cell type
  aggregation_per_celltype <- list()
  # get the unique cell types
  cell_types <- unique(seurat_object@meta.data[[celltype_column]])
  # excluding NA of course
  cell_types <- cell_types[!is.na(cell_types)]
  # check each cell type
  for (cell_type in cell_types) {
    if (verbose) {
      message(paste('calculating for', cell_type))
    }
    # get indices of cells that are of this cell type
    indices_cell_type <- which(seurat_object@meta.data[[celltype_column]] == cell_type)
    # get the count matrix where we have the correct cell type
    count_matrix_full <- GetAssayData(seurat_object[, indices_cell_type], slot = "counts")
    # ignore genes that are never expressed
    count_matrix <-  count_matrix_full[which(rowSums(count_matrix_full) != 0), ]
    
    # also get the metadata for this cell type
    metadata <- seurat_object@meta.data[indices_cell_type, ]
    
    # extract IDs
    ids <- metadata[[participant_column]]
    unique_id_list <- unique(ids)
    
    # create new object to store the counts in
    norm_count_matrix <- count_matrix
    rm(count_matrix)
    gc()
    
    if (verbose) {
      message('doing single-cell normalization')
    }
    # do mean sample-sum normalization
    sample_sum_info = colSums(norm_count_matrix)
    mean_sample_sum = mean(sample_sum_info)
    sample_scale = sample_sum_info / mean_sample_sum
    
    # divide each column by sample_scale
    norm_count_matrix@x <- norm_count_matrix@x / rep.int(sample_scale, diff(norm_count_matrix@p))
    
    if (verbose) {
      message('calculating mean expression')
    }
    # now create the mean expression matrix
    aggregate_norm_count_matrix <- as.data.frame(
      pblapply(
        # go through each ID
        unique_id_list, FUN = function(x){
          # get the sparse means over the cells of a participant
          sparse_Means(norm_count_matrix[, ids == x, drop = FALSE], rowMeans = TRUE)
        }
      )
    )
    # set the colnames to be the participants
    colnames(aggregate_norm_count_matrix) <- unique_id_list
    # and the genes as the rows
    rownames(aggregate_norm_count_matrix) <- rownames(norm_count_matrix)
    
    # save the unfiltered mean expression
    aggregate_norm_count_matrix_unfiltered <- aggregate_norm_count_matrix
    
    # calculate the number of cells per donor
    cell_numbers_per_donor <- data.frame(table(metadata[[participant_column]]))
    # set more descriptive column names
    colnames(cell_numbers_per_donor) <- c('participant', 'number')
    
    # filter on aggregates that are based on a certain number of cells
    samples_with_min_cells <- cell_numbers_per_donor[cell_numbers_per_donor[['number']] >= min_cell_number, 'participant']
    aggregate_norm_count_matrix = aggregate_norm_count_matrix[, which(colnames(aggregate_norm_count_matrix) %in% samples_with_min_cells), drop = F]
    # check if we have any cells left
    if (length(aggregate_norm_count_matrix) != 0 & ncol(aggregate_norm_count_matrix) > 0) {
      # remove genes without any variation
      row_var_info = rowVars(as.matrix(aggregate_norm_count_matrix))
      aggregate_norm_count_matrix = aggregate_norm_count_matrix[which(row_var_info != 0), , drop = F]
      # check if we have any genes left
      if (length(aggregate_norm_count_matrix) != 0 & nrow(aggregate_norm_count_matrix) > 0) {
        if (verbose) {
          message('doing mean expression normalization')
        }
        # do inverse normal transform per gene.
        for (r_i in 1:nrow(aggregate_norm_count_matrix)) {
          aggregate_norm_count_matrix[r_i, ] = qnorm((rank(aggregate_norm_count_matrix[r_i, ], na.last = 'keep')-0.5) / sum(!is.na(aggregate_norm_count_matrix[r_i, ])))
        }
        
        if (verbose) {
          message('performing PCA')
        }
        # Do PCA, and select first x components.
        pc_out = prcomp(t(aggregate_norm_count_matrix))
        # check how many pcs we have
        npcs_present <- ncol(pc_out$x)
        # now subset to that number of pcs if we can
        cov_out <- NULL
        if (npcs <= npcs_present) {
          cov_out = pc_out$x[, 1:npcs]
        }
        else{
          # if we have less pcs we should warn
          warning(paste('requested', as.character(npcs), 'pcs, but only', as.character(npcs_present), 'are present for', cell_type))
          cov_out <- pc_out$x
        }
        if (verbose) {
          message(paste('finished', cell_type))
        }
        # put the results in a list
        aggregate_summary <- list('cell_type' = cell_type, 'expression' = aggregate_norm_count_matrix, 'expression_unfiltered' = aggregate_norm_count_matrix_unfiltered, 'pc' = cov_out)
        # which in turn is put into another list
        aggregation_per_celltype[[cell_type]] <- aggregate_summary
      }
      else {
        message(paste('no genes left after checking for variation between genes for', cell_type, 'skipping cell type'))
      }
    }
    else {
      message(paste('no cells left after min_cells filter for', cell_type, ', skipping cell type'))
    }
  }
  return(aggregation_per_celltype)
}

#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param expression_per_celltype list with celltypes as key, with the filtered expression under the expression key, unfiltered expression under the expression_unfiltered key, and pcs under the pc key
#' @param metadata_per_celltype metadata annation per cell type in a list, where the keys are the cell types
#' @param output_loc where to place the output files
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param rename_samples_to_samplepools rename the sample names 'TEST_81' to the sample+pool 'TEST_81;;190102_lane1'. Required if pools were joined when creating the expression matrices
#' @returns 0 if succesfull
#' write_limix_input(expression_per_celltype, metadata_per_celltype, '/groups/umcg-franke-scrna/tmp02/projects/venema-2022/ongoing/qtl/eqtl/sc-eqtlgen/input/elmentaite_adult_martin_immune/cell_type_safe/')
write_limix_input <- function(expression_per_celltype, metadata_per_celltype, output_loc='./', merge_pcs_into_covariates=F, rename_samples_to_samplepools=T) {
  # we can only do the data that we have expression and metadata for
  cell_types <- intersect(names(expression_per_celltype), names(metadata_per_celltype))
  # create the directory to place the files in if it does not exist yet
  dir.create(output_loc, recursive = T)
  # now check each cell type
  for (cell_type in cell_types) {
    # extract the data
    metadata <- metadata_per_celltype[[cell_type]]
    mean_expression <- expression_per_celltype[[cell_type]][['expression_unfiltered']]
    qtl_expression <- expression_per_celltype[[cell_type]][['expression']]
    pcs <- expression_per_celltype[[cell_type]][['pc']]
    # rename the samples with sample+pool, if not done already
    if (rename_samples_to_samplepools) {
      # subset the expression with entries we have metadata for
      mean_expression <- mean_expression[, colnames(mean_expression) %in% metadata[['IID']]]
      qtl_expression <- qtl_expression[, colnames(qtl_expression) %in% metadata[['IID']]]
      # replace the sample names with sample+lane
      colnames(mean_expression) <- metadata[
        match(colnames(mean_expression), metadata[['IID']]), 'Donor_Pool'
      ]
      colnames(qtl_expression) <- metadata[
        match(colnames(qtl_expression), metadata[['IID']]), 'Donor_Pool'
      ]
      rownames(pcs) <- metadata[
        match(rownames(pcs), metadata[['IID']]), 'Donor_Pool'
      ]
    }
    else{
      # subset the expression with entries we have metadata for
      mean_expression <- mean_expression[, colnames(mean_expression) %in% metadata[['Donor_Pool']]]
      qtl_expression <- qtl_expression[, colnames(qtl_expression) %in% metadata[['Donor_Pool']]]
    }
    # paste together the output names
    qtl_output_loc <- paste(output_loc, '/', cell_type, '.qtlInput.txt', sep = '')
    pcs_output_loc <- paste(output_loc, '/', cell_type, '.qtlInput.Pcs.txt', sep = '')
    exp_output_loc <- paste(output_loc, '/', cell_type, '.Exp.txt', sep = '')
    metadata_output_loc <- paste(output_loc, '/', cell_type, '.covariates.txt', sep = '')
    # write the files
    write.table(qtl_expression, qtl_output_loc, quote = F, sep = '\t', col.names = NA)
    write.table(mean_expression, exp_output_loc, quote = F, sep = '\t', col.names = NA)
    # change X.FFID back to #FID
    colnames(metadata) <- gsub('X\\.FID', '#FID', colnames(metadata))
    # either write the PCs together or separate from the covariates
    if (merge_pcs_into_covariates) {
      # turn into dataframes
      pcs <- data.frame(pcs)
      metadata <- data.frame(metadata)
      # extract samples from PCs
      pc_samples <- rownames(pcs)
      # extract the original column names
      pc_columns <- colnames(pcs)
      # extract the covariate samples
      cov_samples <- as.character(metadata[['Donor_Pool']])
      # extract the original columns
      cov_columns <- colnames(metadata)
      # check which we have in both cases
      joint_samples <- intersect(pc_samples, cov_samples)
      # subset both
      pcs <- pcs[pc_samples, ]
      metadata[as.character(metadata[['Donor_Pool']]) %in% joint_samples, ]
      # add the donor pool to the pcs as a explicit column
      pcs[['Donor_Pool']] <- rownames(pcs)
      # finally join them
      metadata <- merge(metadata, pcs, by = 'Donor_Pool')
      # the reorder back the columns
      metadata <- metadata[, c(cov_columns, pc_columns)]
      # change X.FFID back to #FID
      colnames(metadata) <- gsub('X\\.FID', '#FID', colnames(metadata))
      # and write the result
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = T, row.names = F)
    }
    else {
      write.table(pcs, pcs_output_loc, quote = F, sep = '\t', col.names = NA)
      write.table(metadata, metadata_output_loc, quote = F, sep = '\t', col.names = NA)
    }
  }
  # extract the first metadata
  metadata_first <- metadata_per_celltype[[1]]
  # get the unique combinations of the donor pool and the donor
  smf <-  unique(metadata_first[,c('IID', 'Donor_Pool')])
  colnames(smf) <- c('genotype_id', 'phenotype_id')
  # create the output file
  smf_output_loc <- paste(output_loc, 'smf.txt', sep = '')
  # write the result
  write.table(smf, smf_output_loc, quote=F, sep = '\t', row.names=F)
  return(0)
}

#' get table with the cell numbers for each combination of supplied columns
#' 
#' @param seurat_object_metadata the metadata (from Seurat) to use to create a metadata annotation file
#' @param output_loc where to place the output files
#' @param participant_column the seurat metadata column that denotes the participant
#' @param pool_column the seurat metadata column that denotes the pool that the samples were run in
#' @param condition_column the seurat metadata column that denotes the condition of the sample (inflamed, non-inflamed)
#' @param celltype_column the seurat metadata column that denotes the celltype of the cell
#' @param join_pools whether to join the pools
#' @param min_cell_number the minimal number of cells to need to build a pseudobulk, pseudobulks with less cells are removed
#' @param min_numi the minimal number of UMIs to include a cell for pseudobulk
#' @param npcs the number of PCs to return
#' @param merge_pcs_into_covariates whether to merge the PCs into the covariates file
#' @param verbose print progress or not
#' @returns 0 if succesfull
#' do_limix_input_pipeline(seurat_object, psam, output_loc='./', partipant_column='soup_final_sample_assignment')
do_limix_input_pipeline <- function(seurat_object, 
                                    psam, 
                                    output_loc='./',
                                    participant_column='donor_final', 
                                    pool_column='day', 
                                    condition_column='inflammation_status', 
                                    celltype_column='cell_type_safe', 
                                    join_pools=T,
                                    min_cell_number=5, 
                                    min_numi=200,
                                    npcs=10,
                                    merge_pcs_into_covariates=F, 
                                    verbose=T) {
  if (verbose) {
    message('creating metadata files')
  }
  # create metadata tables
  metadata_per_celltype <- create_metadata_tables(
    seurat_object_metadata = seurat_object@meta.data, 
    psam = psam, 
    participant_column = participant_column, 
    pool_column = pool_column, 
    condition_column = condition_column, 
    celltype_column = celltype_column,
    join_pools = join_pools
  )
  if (verbose) {
    message('creating mean expression matrices')
  }
  # if we joined the pools, we need that as a parameter for the count matrices
  batch_column <- NULL
  # as well as that we don't need to rename the expression matrix columns then
  rename_samples_to_samplepools <- T
  if (!join_pools) {
    batch_column <- pool_column
    rename_samples_to_samplepools <- F
  }
  # create aggregated counts
  expression_per_celltype <- create_aggregated_expression_matrices(
    seurat_object = seurat_object, 
    participant_column = participant_column, 
    celltype_column = celltype_column, 
    batch_column = batch_column,
    min_cell_number = min_cell_number, 
    npcs = npcs, 
    min_numi = min_numi,
    verbose = verbose
  )
  if (verbose) {
    message('writing results')
  }
  # write the results
  write_limix_input(
    expression_per_celltype = expression_per_celltype, 
    metadata_per_celltype = metadata_per_celltype, 
    output_loc = output_loc,
    merge_pcs_into_covariates = merge_pcs_into_covariates, 
    rename_samples_to_samplepools = rename_samples_to_samplepools
  )
  return(0)
}

####################
# Main Code        #
####################

# we will use Seurat version 5
options(Seurat.object.assay.version = 'v5')

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(7777)

# location of the objects
cardio_integrated_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'
m1_v2_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/seurat/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds'
m1_v3_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/seurat/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds'

# age sex file for 1M data
m1_age_sex_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
# age sex file for STEMI data
stemi_age_sex_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/stemi_age_gender_match_wtestid.tsv'
# read the age+sex files
m1_age_sex <- read.table(m1_age_sex_loc, sep = '\t', header = T)
stemi_age_sex <- read.table(stemi_age_sex_loc, sep = '\t', header = T)
# create the m1 psam file
m1_psam <- data.frame(FID = m1_age_sex[['LLD.ID']], IID = m1_age_sex[['LLD.ID']], PAT = 0, MAT = 0, SEX = 0, age = m1_age_sex[['Age']])
m1_psam[m1_age_sex[['Gender']] == 'M', 'SEX'] <- 1
m1_psam[m1_age_sex[['Gender']] == 'F', 'SEX'] <- 2
colnames(m1_psam) <- c('#FID','IID','PAT','MAT','SEX', 'age')
# and the STEMI one
stemi_psam <- data.frame(FID = stemi_age_sex[['ID']], IID = stemi_age_sex[['ID']], PAT = 0, MAT = 0, SEX = 0, age = stemi_age_sex[['age']])
stemi_psam[stemi_age_sex[['gender']] == 'M', 'SEX'] <- 1
stemi_psam[stemi_age_sex[['gender']] == 'F', 'SEX'] <- 2
colnames(stemi_psam) <- c('#FID','IID','PAT','MAT','SEX', 'age')
# and combine them
psam <- rbind(m1_psam, stemi_psam)

# read the stemi object
stemi <- readRDS(cardio_integrated_loc)
# subset to STEMI v2
stemi_v2 <- stemi[, stemi@meta.data[['orig.ident']] == 'stemi_v2']
# clear memory
rm(stemi)
# read the NC2022 data
m1_v2 <- readRDS(m1_v2_loc)
# subset to the unstimulated data
m1_v2_C <- m1_v2[, !is.na(m1_v2@meta.data[['timepoint']]) & m1_v2@meta.data[['timepoint']] == 'UT']
# clear up all the memory again
rm(m1_v2)
# set harmonized column lanes
stemi_v2@meta.data[['sample_condition']] <- stemi_v2@meta.data[['timepoint.final']]
stemi_v2@meta.data[['donor']] <- stemi_v2@meta.data[['assignment.final']]
m1_v2_C[['sample_condition']] <- 'C'
m1_v2_C[['donor']] <- m1_v2_C@meta.data[['assignment']]
# join them together
v2_all <- merge(m1_v2_C, stemi_v2)
# clear memory
rm(stemi_v2)
rm(m1_v2_C)
# get the metadata tables
v2_metadata_tables <- create_metadata_tables(v2_all@meta.data, psam, participant_column='donor', pool_column='lane', condition_column='sample_condition', celltype_column='cell_type_lowerres', join_pools=F)
# get the expression
v2_expression <- create_aggregated_expression_matrices(v2_all, participant_column='donor', celltype_column='cell_type_lowerres', batch_column='lane', min_cell_number=5, min_numi=200, npcs=10, verbose=T)
# and write the results
write_limix_input(expression_per_celltype = v2_expression, metadata_per_celltype = v2_metadata_tables, output_loc = '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v2/', merge_pcs_into_covariates = T, rename_samples_to_samplepools = F)
# remove the object, as we don't need it anymore
rm(v2_all)

# now do everything for v3
stemi <- readRDS(cardio_integrated_loc)
stemi_v3 <- stemi[, stemi@meta.data[['orig.ident']] == 'stemi_v3']
rm(stemi)
m1_v3 <- readRDS(m1_v3_loc)
m1_v3_C <- m1_v3[, !is.na(m1_v3@meta.data[['timepoint']]) & m1_v3@meta.data[['timepoint']] == 'UT']
rm(m1_v3)
stemi_v3@meta.data[['sample_condition']] <- stemi_v3@meta.data[['timepoint.final']]
stemi_v3@meta.data[['donor']] <- stemi_v3@meta.data[['assignment.final']]
m1_v3_C[['sample_condition']] <- 'C'
m1_v3_C[['donor']] <- m1_v3_C@meta.data[['assignment']]
v3_all <- merge(m1_v3_C, stemi_v3)
rm(stemi_v3)
rm(m1_v3_C)
# now all steps at once
do_limix_input_pipeline(v3_all, psam, 
                        output_loc = '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v3/', 
                        participant_column='donor', 
                        pool_column='lane', 
                        condition_column='sample_condition', 
                        celltype_column='cell_type_lowerres', 
                        join_pools=F,
                        min_cell_number=5, 
                        min_numi=200,
                        npcs=10,
                        merge_pcs_into_covariates=T, 
                        verbose=T)
rm(v3_all)