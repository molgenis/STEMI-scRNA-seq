##################################################################################################
## Script for univariate regression analysis on clinical vairables and proportions of monocytes ##
##################################################################################################

library(Seurat)
library(ggplot2)

# get the cell numbers of a cell type and of the total number of cells
get_cell_numbers <- function(metadata, cell_type_to_check, conditions=c('Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', assignment_column='assignment.final' ){
  counts <- NULL
  # check each participant
  for(participant in unique(metadata[[assignment_column]])){
    numbers <- c(participant)
    # check each timepoint
    for(condition in conditions){
      # get the total number of cells in this condition
      total <- nrow(metadata[metadata[[assignment_column]] == participant & metadata[[condition_column]] == condition, ])
      # get the cells of the cell type
      ct <- nrow(metadata[metadata[[assignment_column]] == participant & metadata[[condition_column]] == condition & metadata[[cell_type_column]] == cell_type_to_check, ])
      # add the numbers
      numbers <- c(numbers, total, ct)
    }
    # add to table
    if(is.null(counts)){
      # create the table if non-existant
      counts <- data.frame(t(matrix(numbers)), stringsAsFactors = F)
      each_number <- c('all', cell_type_to_check)
      ct_tp_combinations <- paste(rep(conditions, each = length(each_number)), each_number, sep = "_")
      # put the participant in there, and the columns of each time point wiht 'all' and the cell type we are looking at
      colnames(counts) <- c('participant', ct_tp_combinations)
    }
    else{
      counts <- rbind(counts, numbers)
    }
  }
  return(counts)
}

# load the object
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated_20200820.rds')
# get just the stemi metadata
metadata <- cardio.integrated@meta.data
metadata.stemi <- metadata[metadata$orig.ident == 'stemi_v2' | metadata$orig.ident == 'stemi_v3', ]

# get the cell type numbers from the metadata
cell_type_numbers <- get_cell_numbers(metadata.stemi, 'monocyte')

# read the clinical variable table
clin <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/clinical_parameters/scRNAseq_clinical_data_20200910_final_notextfields.tsv', sep = '\t', header = T, dec = ',')
# Renaming the record ID
clin[clin$record_id == 72,]$record_id <- 73 #apparently already done
# add the 'TEST_' prepend
clin$participant <- paste('TEST_', clin$record_id, sep='')
# merge the two tables
clin_merged <- merge(clin, cell_type_numbers, by='participant')
# the rest of Irene her code expects the record ID to be used, so I will reluctantly change this to be the same as the participant column
clin_merged$record_id <- clin_merged$participant
# and so remove the other column, as it is now the same
clin_merged$participant <- NULL
# Defining v2 and v3 chemistry samples, 1 is V3 and 0 is V2 (Irene's method)
clin_merged$chem <- NA
clin_merged[clin_merged$record_id == "TEST_58" | clin_merged$record_id == "TEST_60" | clin_merged$record_id == "TEST_61"| clin_merged$record_id == "TEST_62" | clin_merged$record_id == "TEST_64" | clin_merged$record_id == "TEST_65" | clin_merged$record_id == "TEST_66" | clin_merged$record_id == "TEST_68" | clin_merged$record_id == "TEST_69" | clin_merged$record_id == "TEST_70" | clin_merged$record_id == "TEST_71" | clin_merged$record_id == "TEST_73" | clin_merged$record_id == "TEST_77" | clin_merged$record_id == "TEST_79" | clin_merged$record_id == "TEST_81" | clin_merged$record_id == "TEST_88", ]$chem <- 1
clin_merged[clin_merged$record_id == "TEST_1" | clin_merged$record_id == "TEST_6" | clin_merged$record_id == "TEST_12" | clin_merged$record_id == "TEST_14" | clin_merged$record_id == "TEST_15" | clin_merged$record_id == "TEST_17" | clin_merged$record_id == "TEST_18" | clin_merged$record_id == "TEST_23" | clin_merged$record_id == "TEST_25" | clin_merged$record_id == "TEST_28" | clin_merged$record_id == "TEST_32" | clin_merged$record_id == "TEST_39" | clin_merged$record_id == "TEST_40" | clin_merged$record_id == "TEST_42" | clin_merged$record_id == "TEST_43" | clin_merged$record_id == "TEST_45" | clin_merged$record_id == "TEST_46"| clin_merged$record_id == "TEST_47"| clin_merged$record_id == "TEST_50" | clin_merged$record_id == "TEST_51"| clin_merged$record_id == "TEST_52"| clin_merged$record_id == "TEST_53"| clin_merged$record_id == "TEST_55"| clin_merged$record_id == "TEST_56", ]$chem <-0

##################################################################################
# Defining clinical variables
##################################################################################
# Independent variables: Pre-PCI (TIMI 0 - TIMI1-3), post-PCI (TIMI0-2 - TIMI3), CRP at baseline, peak CK (log U/L), peak CK-MB (log U/L), peak Troponin T (log U/L), ischemic time (min)
# TIMI flow pre-PCI
clin_merged$preflow <- NA  
clin_merged[!is.na(clin_merged$timi_flow) & clin_merged$timi_flow == 0, ]$preflow <- T
clin_merged[!is.na(clin_merged$timi_flow) & (clin_merged$timi_flow == 1 | clin_merged$timi_flow == 2 | clin_merged$timi_flow == 3), ]$preflow <- F

# Defining TIMI flow pre-PCI as dummy variable
clin_merged$preflow_int <- NA
clin_merged[!is.na(clin_merged$preflow) & clin_merged$preflow == T, ]$preflow_int <- 1
clin_merged[!is.na(clin_merged$preflow) & clin_merged$preflow == F, ]$preflow_int <- 0

# TIMI flow post-PCI
clin_merged$postflow <- NA  
clin_merged[!is.na(clin_merged$timi_flow_postpci) & (clin_merged$timi_flow_postpci == 0 | clin_merged$timi_flow_postpci == 1), ]$postflow <- T
clin_merged[!is.na(clin_merged$timi_flow_postpci) & (clin_merged$timi_flow_postpci == 2 | clin_merged$timi_flow_postpci == 3), ]$postflow <- F

# Defining TIMI flow post-PCI as dummy variable
clin_merged$postflow_int <- NA
clin_merged[!is.na(clin_merged$postflow) & clin_merged$postflow == T, ]$postflow_int <- 1
clin_merged[!is.na(clin_merged$postflow) & clin_merged$postflow == F, ]$postflow_int <- 0

# Defining peak values of CK, CK-MB and Troponin-T
# Peak CK
# Getting lab time when CK was measured
cks <- clin_merged[, c('ck', 'ck_3h', 'ck_6h', 'ck_9h', 'ck_12h', 'ck_24h', 'ck_48h')]
cks[is.na(cks)] <- -1
# Defining new variable
clin_merged$peak_ck <- NA  
# Finding peak value
peak_ck <- apply(cks, 1, max)
# Giving variable a value
clin_merged$peak_ck <- peak_ck
# fix -1 back to NA
clin_merged$peak_ck[clin_merged$peak_ck == -1] <- NA
# Looking at first values of the variable to check if it was succesfull
head(clin_merged[, colnames(clin_merged)[grep("ck", colnames(clin_merged))]])

# Peak CK-MB
# Getting lab time when CK was measured
ck_mbs <- clin_merged[, c('ck_mb', 'ck_mb_3h', 'ck_mb_6h', 'ck_mb_9h', 'ck_mb_12h', 'ck_mb_24h', 'ck_mb_48h')]
ck_mbs[is.na(ck_mbs)] <- -1
# Defining new variable
clin_merged$peak_ckmb <- NA  
# Finding peak value
peak_ckmb <- apply(ck_mbs, 1, max)
# Giving variable a value
clin_merged$peak_ckmb <- peak_ckmb
# fix -1 back to NA
clin_merged$peak_ckmb[clin_merged$peak_ckmb == -1] <- NA
# Looking at first values of the variable to check if it was succesfull
head(clin_merged[, colnames(clin_merged)[grep("ck_mb", colnames(clin_merged))]])

# Peak Troponin T
# Getting lab time when CK was measured
troponints <- clin_merged[, c('troponin_t', 'troponin_t_3h', 'troponin_t_6h', 'troponin_t_9h', 'troponin_t_12h', 'troponin_t_24h', 'troponin_t_48h')]
troponints[is.na(troponints)] <- -1
# Defining new variable
clin_merged$peak_troponint <- NA
# Finding peak value
peak_troponint <- apply(troponints, 1, max)
# Giving variable a value
clin_merged$peak_troponint <- peak_troponint
# fix -1 back to NA
clin_merged$peak_troponint[clin_merged$peak_troponint == -1] <- NA
# Looking at first values of the variable to check if it was succesfull
head(clin_merged[, colnames(clin_merged)[grep("troponin_t", colnames(clin_merged))]])

##################################################################################
# Transforming variables in order to obtain normal distribution
##################################################################################
# Log transforming clinical variables
# CRP
clin_merged$log_crp <- NA
clin_merged$log_crp <- log10(clin_merged$crp)

# Peak-CK
clin_merged$log_peak_ck <- NA  
clin_merged$log_peak_ck <- log10(clin_merged$peak_ck)

# CK-MB
clin_merged$log_peak_ckmb <- NA  
clin_merged$log_peak_ckmb <- log10(clin_merged$peak_ckmb)

# Troponin T
clin_merged$log_peak_troponint <- NA  
clin_merged$log_peak_troponint <- log10(clin_merged$peak_troponint)

# Ischemic time
clin_merged$log_ischemia_time <- NA
clin_merged$log_ischemia_time <- log10(clin_merged$ischemia_time)

# THE REST IS UP TO IRENE...
