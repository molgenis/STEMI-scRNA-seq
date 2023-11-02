#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_pathways_to_excel.R
# Function:
############################################################################################################################

####################
# libraries        #
####################

library(xlsx)


####################
# Functions        #
####################

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
  label_dict[['UT']] <- 'C'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["bulk"]] <- "bulk-like"
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
  return(label_dict)
}


get_pathways_per_celltype_condition <- function(pathway_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), de_method='limma', append='.tsv', pval_column='Adjusted.P.value', pval_cutoff=0.05, direction_column='enrichmentScore', only_positive=F, only_negative=F) {
  # we will save everything in a list first
  output_per_ct_and_comb <- list()
  # we will check each cell type
  for (cell_type in cell_types) {
    # and each condition
    for (condition_1 in conditions) {
      # against each other condition
      for (condition_2 in conditions) {
        # we'll paste the full path together
        path_to_output <- NULL
        if ('limma' == de_method) {
          path_to_output <- paste(pathway_output_loc, '/', cell_type, '_timepoint.final.', condition_1, '_', condition_2, append, sep = '')
        }
        else{
          path_to_output <- paste(pathway_output_loc, '/', cell_type, condition_1, condition_2, append, sep = '')
        }
        # check if the file exists
        if (file.exists(path_to_output)) {
          # read the file
          pathway_output <- read.table(path_to_output, header = T, sep = '\t', comment.char='')
          # filter direction if requested
          if (only_positive) {
            pathway_output <- pathway_output[pathway_output[[direction_column]] > 0, ]
          }
          if (only_negative) {
            pathway_output <- pathway_output[pathway_output[[direction_column]] < 0, ]
          }
          # add the cell type and the condition
          pathway_output <- cbind(
            data.frame(
              condition1=rep(condition_1, times = nrow(pathway_output)),
              condition2=rep(condition_2, times = nrow(pathway_output)),
              cell_type=rep(cell_type, times = nrow(pathway_output))
            ), 
            pathway_output
          )
          # filter if possible
          if (!is.null(pval_column) & !is.null(pval_cutoff)) {
            pathway_output <- pathway_output[pathway_output[[pval_column]] < pval_cutoff, ]
          }
          # add to the list
          output_per_ct_and_comb[[paste(condition_1, condition_2, cell_type, sep = '')]] <- pathway_output
        }
      }
    }
  }
  # merge everything
  output_per_ct_and_comb_df <- do.call('rbind', output_per_ct_and_comb)
  return(output_per_ct_and_comb_df)
}


pathways_to_excel <- function(pathway_table, excel_output_loc, cell_type_column='cell_type', direction_column='direction', condition_1_column='condition1', condition_2_column='condition2' ,use_label_dict=T, remove_columns=c()) {
  # create a new workbook
  wb = createWorkbook()
  # get the cell types
  cell_types <- unique(pathway_table[[cell_type_column]])
  # check each cell type
  for (celltype in cell_types) {
    # get the table for that cell type
    pathway_table_celltype <- pathway_table[!is.na(pathway_table[[cell_type_column]]) & pathway_table[[cell_type_column]] == celltype, ]
    # create the name for the sheet
    sheet_name <- celltype
    # replace with a better label if requested
    if (use_label_dict) {
      sheet_name <- label_dict()[[celltype]]
    }
    # check specifically the up and down
    directions <- unique(pathway_table_celltype[[direction_column]])
    for (direction in directions) {
      # subset to that direction
      pathway_table_celltype_direction <- pathway_table_celltype[!is.na(pathway_table_celltype[[direction_column]]) & pathway_table_celltype[[direction_column]] == direction, ]
      # paste together the tab name
      tab_name <- paste(sheet_name, direction)
      # create the sheet
      sheet = createSheet(wb, tab_name)
      # again, replace labels if requested
      if (use_label_dict) {
        pathway_table_celltype_direction[[cell_type_column]] <- as.vector(unlist(label_dict()[pathway_table_celltype_direction[[cell_type_column]]]))
        pathway_table_celltype_direction[[condition_1_column]] <- as.vector(unlist(label_dict()[pathway_table_celltype_direction[[condition_1_column]]]))
        pathway_table_celltype_direction[[condition_2_column]] <- as.vector(unlist(label_dict()[pathway_table_celltype_direction[[condition_2_column]]]))
      }
      if (!is.null(remove_columns)) {
        pathway_table_celltype_direction <- pathway_table_celltype_direction[, setdiff(colnames(pathway_table_celltype_direction), remove_columns)]
      }
      # add dataframe to sheet
      addDataFrame(pathway_table_celltype_direction, sheet = sheet, startColumn=1, row.names=FALSE)
    }
  }
  # write the result
  saveWorkbook(wb, excel_output_loc)
  return(0)
}

####################
# Main Code        #
####################

# get the location of the pathway output
pathway_output_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/pathway_enrichment/enrichr_reactome_2023/limma_dream/stemi_both_paired_lores_lfc01minpct01ncountrna_20210301/sct/chem_assignment_age_gender/'

# get all output
limma_both_enrichr_negative <- get_pathways_per_celltype_condition(pathway_output_loc, pval_cutoff = 0.05, pval_column = 'Old.Adjusted.P.value', append = '_negative.tsv')
limma_both_enrichr_positive <- get_pathways_per_celltype_condition(pathway_output_loc, pval_cutoff = 0.05, pval_column = 'Old.Adjusted.P.value', append = '_positive.tsv')

# add direction as explicit column
limma_both_enrichr_negative <- cbind(data.frame(direction = rep('down', times = nrow(limma_both_enrichr_negative))), limma_both_enrichr_negative)
limma_both_enrichr_positive <- cbind(data.frame(direction = rep('up', times = nrow(limma_both_enrichr_positive))), limma_both_enrichr_positive)

# and add together
limma_both_enrichr <- rbind(limma_both_enrichr_positive, limma_both_enrichr_negative)

# where we store the Excel
limma_both_enrichr_excel_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/pathway_enrichment/enrichr_reactome_2023/stemi_limma_pathways_reactome_2023.xlsx'

# create the Excel
pathways_to_excel(limma_both_enrichr, limma_both_enrichr_excel_loc)
