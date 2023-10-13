#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_meta_analys_de.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################

library(metaRNASeq)

####################
# Functions        #
####################


add_ncell_ndonor <- function(mast_output_loc, mast_added_output_loc, metadata_loc='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/cardio.integrated.20210301.metadata.tsv', cell_type_column='cell_type_lowerres', timepoint_column='timepoint.final', participant_column='assignment.final', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), timepoints=c('UT', 'Baseline', 't24h', 't8w')) {
  # read the metadata file
  metadata <- read.table(metadata_loc, header = T, sep = '\t')
  # create a table with cell numbers
  cell_numbers <- data.frame(table(metadata[, c(cell_type_column, timepoint_column, participant_column)]))
  # remove empty entries
  cell_numbers <- cell_numbers[cell_numbers[['Freq']] > 0, ]
  # save data in a list
  metadata_numbers <- list()
  # check each timepoint
  for (timepoint in unique(cell_numbers[[timepoint_column]])) {
    # add timepoint to the list
    metadata_numbers[[timepoint]] <- list()
    # subset to the timepoint
    cell_numbers_timepoint <- cell_numbers[cell_numbers[[timepoint_column]] == timepoint, ]
    # check each cell type
    for (cell_type in unique(cell_numbers_timepoint[[cell_type_column]])) {
      # subset to the cell type
      cell_numbers_timepoint_celltype <- cell_numbers_timepoint[cell_numbers_timepoint[[cell_type_column]] == cell_type, ]
      # get the number of cells
      ncell <- sum(cell_numbers_timepoint_celltype[['Freq']])
      # get the number of donors
      ndonor <- length(unique(cell_numbers_timepoint_celltype[[participant_column]]))
      # add this to the list
      metadata_numbers[[timepoint]][[cell_type]] <- list('ncell' = ncell, 'ndonor' = ndonor)
    }
  }
  # now let's check our DE output
  for (cell_type in cell_types) {
    # and each timepoint
    for (timepoint1 in timepoints) {
      # against each timepoint
      for(timepoint2 in timepoints) {
        # paste together what the path should be
        full_de_path <- paste(mast_output_loc, '/', cell_type, timepoint1, timepoint2, '.tsv', sep = '')
        # check if that file exists
        if (file.exists(full_de_path)) {
          # read the file
          de_output <- read.table(full_de_path, header = T, row.names = 1, sep = '\t')
          # extract the genes
          genes <- rownames(de_output)
          # get the number of cells in timeopint 1
          ncell_1 <- metadata_numbers[[timepoint1]][[cell_type]][['ncell']]
          # and timepoint 2
          ncell_2 <- metadata_numbers[[timepoint2]][[cell_type]][['ncell']]
          # and donors
          ndonor_1 <- metadata_numbers[[timepoint1]][[cell_type]][['ndonor']]
          ndonor_2 <- metadata_numbers[[timepoint2]][[cell_type]][['ndonor']]
          # and comparison
          comparison <- paste(timepoint1, timepoint2, sep = '-')
          # make the stat df
          stat_info <- data.frame(
            gene = genes,
            combination = rep(comparison, times = length(genes)),
            ndonor_t1 = rep(ndonor_1, times = length(genes)),
            cell_total_t1 <- rep(ncell_1, times = length(genes)),
            ndonor_t2 = rep(ndonor_2, times = length(genes)),
            cell_total_t2 <- rep(ncell_2, times = length(genes))
          )
          # set the column names as expected
          colnames(stat_info) <- c(
            'gene',
            'combination',
            paste('ndonor_', timepoint1, sep = ''),
            paste('cell_total_', timepoint1, sep = ''),
            paste('ndonor_', timepoint2, sep = ''),
            paste('cell_total_', timepoint2, sep = '')
          )
          # merge with original output
          de_output <- cbind(stat_info, de_output)
          # set the rownames back
          rownames(de_output) <- genes
          # create the output path
          new_de_path <- paste(mast_added_output_loc, '/', cell_type, timepoint1, timepoint2, '.tsv', sep = '')
          # write output
          write.table(de_output, new_de_path, row.names = T, col.names = T, sep = '\t')
        }
      }
    }
  }
  return(0)
}


write_meta_de <- function(de_output_loc_prepend, de_output_loc_append, de_meta_output_loc_prepend, sample_number_regex='ndonor_*', pval_column='P.Value', fc_column='logFC', file_regex='*.tsv'){
  # go through the conditions
  files <- list.files(paste(de_output_loc_prepend, '2', de_output_loc_append, sep = ''), pattern = file_regex)
  for(file in files){
    # get the de output
    de_loc_v2 <- paste(de_output_loc_prepend, '2', de_output_loc_append, file, sep = '')
    # v3 should be the same, but then v3
    de_loc_v3 <- paste(de_output_loc_prepend, '3', de_output_loc_append, gsub('v2', 'v3', file), sep = '')
    tryCatch({
      # read the de output
      de_v2 <- read.table(de_loc_v2, header=T)
      de_v3 <- read.table(de_loc_v3, header=T)
      # get the genes that are in both
      genes_both <- intersect(rownames(de_v2), rownames(de_v3))
      # select only those genes
      de_v2 <- de_v2[genes_both, ]
      de_v3 <- de_v3[genes_both, ]
      # morph P val to minimum
      if(nrow(de_v2[de_v2[[pval_column]] == 0, ]) > 0){
        de_v2[de_v2[[pval_column]] == 0, 'pval_column'] <- .Machine$double.xmin
      }
      if(nrow(de_v3[de_v3[[pval_column]] == 0, ]) > 0){
        de_v3[de_v3[[pval_column]] == 0, 'pval_column'] <- .Machine$double.xmin
      }
      # add the gene name also in a column
      de_v2$gene <- rownames(de_v2)
      de_v3$gene <- rownames(de_v3)
      # get the P values
      rawpval <- list('v2' = de_v2[[pval_column]], 'v3' = de_v3[[pval_column]])
      
      # get the total number of samples
      cols_samples_v2 <- colnames(de_v2)[grep(sample_number_regex, colnames(de_v2))]
      cols_samples_v3 <- colnames(de_v3)[grep(sample_number_regex, colnames(de_v3))]
      # now I have the columns, I just need the first entries
      nsamples <- list(
        'v2' = sum(de_v2[1, cols_samples_v2]),
        'v3' = sum(de_v2[2, cols_samples_v2])
      )
      
      # do fisher combining
      fisher_output <-fishercomb(rawpval)
      
      # and inverse norm
      invnorm_output <- invnorm(rawpval, as.vector(unlist(nsamples)))
      
      # now add everything together
      colnames(de_v2) <- paste('v2', colnames(de_v2), sep = '.')
      colnames(de_v3) <- paste('v3', colnames(de_v3), sep = '.')
      de_meta <- cbind(de_v2, de_v3)
      de_meta[['fisher.TestStatistic']] <- fisher_output[['TestStatistic']]
      de_meta[['fisher.pval']] <- fisher_output[['adjpval']]
      de_meta[['invnorm.TestStatistic']] <- invnorm_output[['TestStatistic']]
      de_meta[['invnorm.pval']] <- invnorm_output[['adjpval']]
      # check the sign
      de_meta[['same_sign']] <- apply(de_meta, 1, function(x) {
        sign(as.numeric(x[paste('v2', fc_column, sep = '.')])) == sign(as.numeric(x[paste('v3', fc_column, sep = '.')]))
        }
      )
      # add a final p
      de_meta[['invnorm.TestStatistic.signed']] <- de_meta[['invnorm.TestStatistic']]
      de_meta[de_meta[['same_sign']] == F, 'invnorm.TestStatistic.signed'] <- NA
      # add meta fc
      de_meta[['meta_fc']] <- apply(de_meta, 1, function(x) {
        mean(as.numeric(x[paste('v2', fc_column, sep = '.')]), as.numeric(x[paste('v3', fc_column, sep = '.')]))
        }
      )
      # write the result
      output_loc <- paste(de_meta_output_loc_prepend, file, sep = '')
      write.table(de_meta, output_loc, sep = '\t')
    }, error = function(e) {
      print(paste('error in', file, 'due to', e))
    })
  }
}

####################
# Main Code        #
####################

# paste directories together
de_output_prepend <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/limma_dream/stemi_v'
de_output_append <- '_paired_lores_lfc01minpct01ncountrna_20210301/sct/assignment_age_gender/'
de_meta_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/limma_dream/stemi_meta_paired_lores_lfc01minpct01ncountrna_20210301/sct/chem_assignment_age_gender/'

# write meta output
dir.create(de_meta_loc, recursive = T)
write_meta_de(de_output_prepend, de_output_append, de_meta_loc)


# add some info to the MAST output
add_ncell_ndonor('/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_v2_paired_lores_lfc01minpct01ncountrna_20210301/rna/', '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_v2_paired_lores_lfc01minpct01ncountrna_20210301_wstats/rna/')
add_ncell_ndonor('/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_v3_paired_lores_lfc01minpct01ncountrna_20210301/rna/', '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_v3_paired_lores_lfc01minpct01ncountrna_20210301_wstats/rna/')
# now do the same for the MAST output
de_mast_output_prepend <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_v'
de_mast_output_append <- '_paired_lores_lfc01minpct01ncountrna_20210301_wstats/rna/'
de_mast_meta_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_metasamesign_paired_lores_lfc01minpct01ncountrna_20210301_wstats_byncell/rna/'
dir.create(de_mast_meta_loc, recursive = T)
write_meta_de(de_mast_output_prepend, de_mast_output_append, de_mast_meta_loc, sample_number_regex = 'cell_total_*', pval_column = 'p_val', fc_column = 'avg_logFC')
