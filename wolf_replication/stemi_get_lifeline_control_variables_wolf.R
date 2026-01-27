#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_get_lifeline_control_variables_wolf.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################


####################
# Functions        #
####################


ll_to_pseudo_ext <- function(ll_to_pseudo, pseudo_int_to_ext){
  # match one onto the other
  ll_to_pseudo$V3 <- pseudo_int_to_ext[match(ll_to_pseudo$V2, pseudo_int_to_ext$PROJECT_PSEUDO_ID), 'PSEUDOIDEXT']
  # set column names
  colnames(ll_to_pseudo) <- c('ll', 'psext', 'psint')
  return(ll_to_pseudo)
}


pseudo_int_to_ugli <- function(ll_to_pseudo, pseudo_int_to_ugli){
  ll_to_pseudo$ugli <- pseudo_int_to_ugli[match(ll_to_pseudo$psint, pseudo_int_to_ugli$PSEUDOIDEXT), 'UGLI_ID']
  return(ll_to_pseudo)
}


pseudo_ext_to_cytoid <- function(table_w_pseudoint, pseudo_int_to_cyto){
  table_w_pseudoint$cyto <- pseudo_int_to_cyto[match(table_w_pseudoint$psint, pseudo_int_to_cyto$PSEUDOIDEXT), 'cytosnp_ID']
  return(table_w_pseudoint)
}


add_metadata <- function(table_w_pseudoint, variables_per_file){
  
}

calculate_stats <- function(pseudo_ext_ids, table_loc, table_column, stats=c('mean', 'median'), id_column='PSEUDOIDEXT'){
  # read the table
  meta_table <- read.table(table_loc, sep = '\t', header = T)
  # grab the interested column
  variable_for_stats <- meta_table[match(pseudo_ext_ids, meta_table[[id_column]]), table_column]
  # set the outputs
  stat_results <- list()
  # go through the stats
  for(stat in stats){
    # if mean
    if(stat == 'mean'){
      stat_results[[stat]] <- mean(variable_for_stats)
    }
    else if(stat == 'median'){
      stat_results[[stat]] <- median(variable_for_stats)
    }
    else if(stat == 'pct'){
      
      if(length(unique(variable_for_stats)) > 2){
        # we only do one or the other
        print('more than two variables, skipping')
        print((variable_for_stats))
      }
      else if(length(unique(variable_for_stats)) < 1){
        print('less than one variable, skipping')
        print((variable_for_stats))
      }
      else{
        # number of unique values
        unique_values <- unique(variable_for_stats)
        # tatal number of values
        total_values <- length(variable_for_stats)
        # amount of values, option 1
        values_1 <- sum(variable_for_stats == unique_values[1])
        # calculate pct
        pct <- values_1/total_values
        stat_results[[paste(stat,'_',unique_values[1],sep = '')]] <- pct
      }
    }
  }
  return(stat_results)
}

stats_to_table <- function(pseudo_ext_ids, metadata_vars, id_column='PSEUDOIDEXT'){
  # init the table we want
  metadata_result <- NULL
  # check each variable
  for(var in names(metadata_vars)){
    # get that entry
    entry <- metadata_vars[[var]]
    # get specific column
    column <- entry[['column']]
    # get specific file
    table_loc <- entry[['table']]
    # get interested stats
    stats <- entry[['stats']]
    # grab the result
    stats <- calculate_stats(pseudo_ext_ids, table_loc, column, stats, id_column)
    # check each calculated stat
    for(stat_name in names(stats)){
      # turn into a row
      row <- data.frame(variable=c(var), stat=c(stat_name), value=c(as.character(stats[[stat_name]])))
      # add to the rest
      if(is.null(metadata_result)){
        metadata_result <- row
      }
      else{
        metadata_result <- rbind(metadata_result, row)
      }
    }
  }
  return(metadata_result)
}


test_expression_to_metadata <- function(expression_table, metadata_table, genes, meta_variables, exp_to_meta_participant_mapping){
  # these genes we will use
  genes_to_use <- intersect(rownames(expression_table), genes)
  # this metadata we will use
  metadata_variables <- intersect(colnames(metadata_table), meta_variables)
  # create a result table
  results <- data.frame(NA, nrow=length(genes_to_use)*length(metadata_variables), ncol=5, dimnames = list(NA, c('gene', 'variable', 'n', 'p', 'method')))
  # get the participants in the meta data
  meta_participants <- metadata_table[['PSEUDOIDEXT']]
  # get the expression participants
  exp_participants <- colnames(expression_table)
  # subset to what is in both of the tables
  mapping_subset <- exp_to_meta_participant_mapping[exp_to_meta_participant_mapping$ll %in% exp_participants &
                                                      exp_to_meta_participant_mapping$psint %in% meta_participants, ]
  
  # row index
  i <- 1
  # check each gene
  for(gene in genes_to_use){
    # check each metadata variable
    for(variable in metadata_variables){
      # grab the expression values
      expression <- as.vector(unlist(expression_table[gene, mapping_subset$ll]))
      
      # grab the metadata 
      metadata_subtable <- metadata_table[match(mapping_subset$psint, metadata_table$PSEUDOIDEXT), ]
      metadata <- metadata_subtable[[variable]]
      # check how many entries
      n <- length(metadata)
      # init entry
      p <- NA
      method <- NA
      # check which analysis to do
      if(is.numeric(metadata) & length(unique(metadata)) > 2){
        # will do correlation
        try({
          p <- cor.test(metadata, expression, method = 'spearman')$p.value
          method = 'spearman'
        })
      }
      else if(is.factor(metadata) | is.character(metadata) | (is.numeric(metadata) & length(unique(metadata)) == 2)){
        # test for normal distribution
        normal <- shapiro.test(expression)
        # do different test depending on normality
        if(normal$p.value < 0.05){
          # non-normal
          wilcox <- wilcox.test(expression, metadata)
          p <- wilcox$p.value
          method <- 'wilcoxon'
        }
        else if(normal$p.value >= 0.05){
          # normal
          ttest <- t.test(expression, metadata)
          p <- ttest$p.value
          method <- 'ttest'
        }
      }
      else{
        print('not numeric, character or factor, or only all identical values')
      }
      # add results
      row <- c(gene, variable, method, p, n)
      results[i, ] <- row
      # increase the index
      i <- i + 1
    }
  }
  return(results)
}


get_relevant_metadata <- function(variables_per_table, interested_pseudo_ids, sep = ','){
  # create the result table
  result_table <- NULL
  # check each variable
  for(table in names(variables_per_table)){
    # get the table
    table_path <- variables_per_table[[table]][['location']]
    # get the interesting columns
    interested_columns <- variables_per_table[[table]][['variables']]
    try({
      # read the table
      table_meta <- read.table(table_path, header = T, stringsAsFactors = F, sep = sep)
      # subset to the participants
      table_meta <- table_meta[match(interested_pseudo_ids, table_meta[['project_pseudo_id']]), ,drop = F]
      # get the overlapping columns
      columns_overlapping <- intersect(colnames(table_meta), interested_columns)
      if(is.null(columns_overlapping) | length(columns_overlapping) == 0) {
        print('no overlap')
      } 
      else {
        # report on what is missing
        missing_columns <- setdiff(interested_columns, colnames(table_meta))
        if(length(missing_columns) > 0){
          print('missing columns:')
          print(paste(missing_columns, collapse = ','))
        }
        # grab the interested columns
        table_interested <- table_meta[, columns_overlapping, drop = F]
        # get the enumerations
        table_enums_path <- variables_per_table[[table]][['enumerations']]
        table_enums <- read.table(table_enums_path, sep = ',', header = T)
        print(head(table))
        # check each column
        for (column in columns_overlapping) {
          # only replace if there is an enumeration
          if (column %in% table_enums[['variable_name']]) {
            # subset to that variable
            enums_column <- table_enums[table_enums[['variable_name']] == column, c('enumeration_code', 'enumeration_en')]
            # check if the enumeration is in the enums
            enums_present <- enums_column[['enumeration_code']]
            # now replace the code with the text value
            print(table_interested[!is.na(table_interested[[column]]) &
                                     table_interested[[column]] %in% enums_present, 
            ])
            table_interested[!is.na(table_interested[[column]]) &
                               table_interested[[column]] %in% enums_present, 
                             column] <- enums_column[match(table_interested[!is.na(table_interested[[column]]) &
                                                                              table_interested[[column]] %in% enums_present, 
                                                                            column], enums_column[['enumeration_code']]), 'enumeration_en']
          }
          else{
          }
        }
        # add the name of the table to the columns
        colnames(table_interested) <- paste(table, '.', colnames(table_interested), sep = '')
        # add to result table
        if(is.null(result_table)){
          result_table <- table_interested
        }
        else{
          result_table <- cbind(result_table, table_interested)
        }
      }
    })
  }
  return(result_table)
}


get_relevant_descriptions <- function(variables_per_table){
  # create the result table
  result_table <- NULL
  # check each variable
  for(variables in variables_per_table){
    # get the table
    table_path <- variables[['tbl_loc']]
    # get the interesting columns
    interested_columns <- variables[['columns']]
    try({
      # read the table
      table_meta <- read.table(table_path, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
      # grab the interested columns
      table_interested <- table_meta[interested_columns, 'Study.subject.identifier', drop = F]
      # add to result table
      if(is.null(result_table)){
        result_table <- table_interested
      }
      else{
        result_table <- rbind(result_table, table_interested)
      }
    })
    
  }
  return(result_table)
}



# add the stimulation tags, so UT, x hours after stim y, etc.
check_missingnes <- function(stim_mapping_loc, ids_available){
  # grab the timepoints
  tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
  # we know now how many entries we will have
  presence_per_lane <- matrix(, ncol=2, nrow=nrow(tps), dimnames = list(rownames(tps), c('n_total', 'n_present')))
  # check the timepoints
  for(lane in rownames(tps)){
    # initialize counter
    nr_present_lane <- 0
    nr_total_lane <- 0
    # check the lanes
    for(tp in colnames(tps)){
      # only apply if there are actually rows with lane in the object
      participants.as.string <- tps[lane,tp]
      # split the participant line by comma to get the participants for the timepoint
      participants.this.tp <- strsplit(participants.as.string, ",")[[1]]
      # only if present we can continue
      if(!is.na(participants.this.tp)){
        # check how many are in the ids vector
        present <- sum(participants.this.tp %in% as.character(ids_available))
        # add to the counter for this lane
        nr_present_lane <- nr_present_lane + present
        # and the overall counter
        nr_total_lane <- nr_total_lane + length(participants.this.tp)
      }
    }
    presence_per_lane[lane, 'n_total'] <- nr_total_lane
    presence_per_lane[lane, 'n_present'] <- nr_present_lane
  }
  return(presence_per_lane)
}


check_missingnes_wmatching <- function(stim_mapping_loc, ids_available, mapping_table=NA, from_column=NA, to_column=NA){
  # by default we will just check the given IDs
  ids_to_check <- ids_available
  # but we will convert if necessary
  if(!is.na(mapping_table) & !is.na(from_column) &!is.na(to_column)){
    # grab the timepoints
    tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
    # set the ids we want
    ids_to_check <- mapping_table[match(ids_to_check, mapping_table[[from_column]]), 'to_column']
  }
  # now do the actual analysis
  presence_per_lane <- check_tags(stim_mapping_loc, ids_available)
  return(presence_per_lane)
}

####################
# Main Code        #
####################

####################
# Fetch all IDs    #
####################

# location of the LL to pseudo
full_id_table_controls_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/ongoing/metadata/stemi_lld_psints.tsv'

# location of the questionaire data
lifelines_tables_loc <- '/groups/umcg-lifelines/rsc02/releases/pheno_lifelines/v2/results/'
lifelines_tables_1a_q_1_results_loc <- paste(lifelines_tables_loc, '1a_q_1_results.csv', sep = '')
lifelines_tables_1a_q_2_results_loc <- paste(lifelines_tables_loc, '1a_q_2_results.csv', sep = '')
lifelines_tables_1a_q_youth_results_loc <- paste(lifelines_tables_loc, '1a_q_youth_results.csv', sep = '')
lifelines_tables_1a_v_1_results_loc <- paste(lifelines_tables_loc, '1a_v_1_results.csv', sep = '')
lifelines_tables_1a_v_2_results_loc <- paste(lifelines_tables_loc, '1a_v_2_results.csv', sep = '')
lifelines_tables_1b_q_1_results_loc <- paste(lifelines_tables_loc, '1b_q_1_results.csv', sep = '')
lifelines_tables_1c_q_1_results_loc <- paste(lifelines_tables_loc, '1c_q_1_results.csv', sep = '')

lifelines_tables_2a_q_1_results_loc <- paste(lifelines_tables_loc, '2a_q_1_results.csv', sep = '')
lifelines_tables_2a_q_2_results_loc <- paste(lifelines_tables_loc, '2a_q_2_results.csv', sep = '')
lifelines_tables_2a_q_youth_results_loc <- paste(lifelines_tables_loc, '2a_q_youth_results.csv', sep = '')
lifelines_tables_2a_v_1_results_loc <- paste(lifelines_tables_loc, '2a_v_1_results.csv', sep = '')
lifelines_tables_2a_v_2_results_loc <- paste(lifelines_tables_loc, '2a_v_2_results.csv', sep = '')
lifelines_tables_2b_q_1_results_loc <- paste(lifelines_tables_loc, '2b_q_1_results.csv', sep = '')

lifelines_tables_3a_q_1_results_loc <- paste(lifelines_tables_loc, '3a_q_1_results.csv', sep = '')
lifelines_tables_3a_q_2_results_loc <- paste(lifelines_tables_loc, '3a_q_2_results.csv', sep = '')
lifelines_tables_3a_q_youth_results_loc <- paste(lifelines_tables_loc, '3a_q_youth_results.csv', sep = '')
lifelines_tables_3a_v_1_results_loc <- paste(lifelines_tables_loc, '3a_v_1_results.csv', sep = '')
lifelines_tables_3b_q_1_results_loc <- paste(lifelines_tables_loc, '3b_q_1_results.csv', sep = '')

# location of the questionaire enumerations
lifelines_tables_loc <- '/groups/umcg-lifelines/rsc02/releases/pheno_lifelines/v2/enumerations/'
lifelines_tables_1a_q_1_enumerations_loc <- paste(lifelines_tables_loc, '1a_q_1_enumerations.csv', sep = '')
lifelines_tables_1a_q_2_enumerations_loc <- paste(lifelines_tables_loc, '1a_q_2_enumerations.csv', sep = '')
lifelines_tables_1a_q_youth_enumerations_loc <- paste(lifelines_tables_loc, '1a_q_youth_enumerations.csv', sep = '')
lifelines_tables_1a_v_1_enumerations_loc <- paste(lifelines_tables_loc, '1a_v_1_enumerations.csv', sep = '')
lifelines_tables_1a_v_2_enumerations_loc <- paste(lifelines_tables_loc, '1a_v_2_enumerations.csv', sep = '')
lifelines_tables_1b_q_1_enumerations_loc <- paste(lifelines_tables_loc, '1b_q_1_enumerations.csv', sep = '')
lifelines_tables_1c_q_1_enumerations_loc <- paste(lifelines_tables_loc, '1c_q_1_enumerations.csv', sep = '')

lifelines_tables_2a_q_1_enumerations_loc <- paste(lifelines_tables_loc, '2a_q_1_enumerations.csv', sep = '')
lifelines_tables_2a_q_2_enumerations_loc <- paste(lifelines_tables_loc, '2a_q_2_enumerations.csv', sep = '')
lifelines_tables_2a_q_youth_enumerations_loc <- paste(lifelines_tables_loc, '2a_q_youth_enumerations.csv', sep = '')
lifelines_tables_2a_v_1_enumerations_loc <- paste(lifelines_tables_loc, '2a_v_1_enumerations.csv', sep = '')
lifelines_tables_2a_v_2_enumerations_loc <- paste(lifelines_tables_loc, '2a_v_2_enumerations.csv', sep = '')
lifelines_tables_2b_q_1_enumerations_loc <- paste(lifelines_tables_loc, '2b_q_1_enumerations.csv', sep = '')

lifelines_tables_3a_q_1_enumerations_loc <- paste(lifelines_tables_loc, '3a_q_1_enumerations.csv', sep = '')
lifelines_tables_3a_q_2_enumerations_loc <- paste(lifelines_tables_loc, '3a_q_2_enumerations.csv', sep = '')
lifelines_tables_3a_q_youth_enumerations_loc <- paste(lifelines_tables_loc, '3a_q_youth_enumerations.csv', sep = '')
lifelines_tables_3a_v_1_enumerations_loc <- paste(lifelines_tables_loc, '3a_v_1_enumerations.csv', sep = '')
lifelines_tables_3b_q_1_enumerations_loc <- paste(lifelines_tables_loc, '3b_q_1_enumerations.csv', sep = '')


# add the cyto ID identifiers as well
full_id_table_X1 <- read.table(full_id_table_controls_loc, header = T, sep = '\t')

####################
# get the ll data  #
####################

# the table that has the names of the variable I want
variable_mapping_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/stemi_lifelines_variables_codes_wolf.tsv'
variable_mapping_variables <- read.table(variable_mapping_loc, sep = '\t', header = F)
# I don't want to look up what is available where, so I'll just try to get it wherever it is available
variable_mapping_list <- c()
for (tbl in c('1aq1', '1aq2', '1aqyouth', '1av1', '1av2', '1bq1', '1cq1', 
              '2aq1', '2aq2', '2aqyouth', '2av1', '2av2', '2bq1', 
              '3aq1', '3aq2', '3aqyouth', '3av1', '3bq1')) {
  # make a table with these variables
  variable_mapping_tbl <- data.frame(
    'V1' = variable_mapping_variables$V1, 
    'V2' = variable_mapping_variables$V1, 
    'V3' = variable_mapping_variables$V1, 
    'V4' = rep(tbl, times = length(variable_mapping_variables$V1))
  )
  # put in list
  variable_mapping_list[[tbl]] <- variable_mapping_tbl
}
# merge all
variable_mapping <- do.call('rbind', variable_mapping_list)
# get the variables for a1
variables_1aq1 <- variable_mapping[variable_mapping[['V4']] == '1aq1', 'V3']
variables_1aq2 <- variable_mapping[variable_mapping[['V4']] == '1aq2', 'V3']
variables_1aqyouth <- variable_mapping[variable_mapping[['V4']] == '1aqyouth', 'V3']
variables_1av1 <- variable_mapping[variable_mapping[['V4']] == '1av1', 'V3']
variables_1av2 <- variable_mapping[variable_mapping[['V4']] == '1av2', 'V3']
variables_1bq1 <- variable_mapping[variable_mapping[['V4']] == '1bq1', 'V3']
variables_1cq1 <- variable_mapping[variable_mapping[['V4']] == '1cq1', 'V3']

variables_2aq1 <- variable_mapping[variable_mapping[['V4']] == '2aq1', 'V3']
variables_2aq2 <- variable_mapping[variable_mapping[['V4']] == '2aq2', 'V3']
variables_2aqyouth <- variable_mapping[variable_mapping[['V4']] == '2aqyouth', 'V3']
variables_2av1 <- variable_mapping[variable_mapping[['V4']] == '2av1', 'V3']
variables_2av2 <- variable_mapping[variable_mapping[['V4']] == '2av2', 'V3']
variables_2bq1 <- variable_mapping[variable_mapping[['V4']] == '2bq1', 'V3']

variables_3aq1 <- variable_mapping[variable_mapping[['V4']] == '3aq1', 'V3']
variables_3aq2 <- variable_mapping[variable_mapping[['V4']] == '3aq2', 'V3']
variables_3aqyouth <- variable_mapping[variable_mapping[['V4']] == '3aqyouth', 'V3']
variables_3av1 <- variable_mapping[variable_mapping[['V4']] == '3av1', 'V3']
variables_3bq1 <- variable_mapping[variable_mapping[['V4']] == '3bq1', 'V3']


# have the variables per result file from lifelines
table_and_columns <- list()
table_and_columns[['1aq1']] <- list('location' = lifelines_tables_1a_q_1_results_loc,  'variables' =  variables_1aq1, 'enumerations' = lifelines_tables_1a_q_1_enumerations_loc)
table_and_columns[['1aq2']] <- list('location' = lifelines_tables_1a_q_2_results_loc,  'variables' =  variables_1aq2, 'enumerations' = lifelines_tables_1a_q_2_enumerations_loc)
table_and_columns[['1aqyouth']] <- list('location' = lifelines_tables_1a_q_youth_results_loc,  'variables' =  variables_1aqyouth, 'enumerations' = lifelines_tables_1a_q_youth_enumerations_loc)
table_and_columns[['1av1']] <- list('location' = lifelines_tables_1a_v_1_results_loc,  'variables' =  variables_1av1, 'enumerations' = lifelines_tables_1a_v_1_enumerations_loc)
table_and_columns[['1av2']] <- list('location' = lifelines_tables_1a_v_2_results_loc,  'variables' =  variables_1av2, 'enumerations' = lifelines_tables_1a_v_1_enumerations_loc)
table_and_columns[['1bq1']] <- list('location' = lifelines_tables_1b_q_1_results_loc,  'variables' =  variables_1bq1, 'enumerations' = lifelines_tables_1b_q_1_enumerations_loc)
table_and_columns[['1cq1']] <- list('location' = lifelines_tables_1b_q_1_results_loc,  'variables' =  variables_1cq1, 'enumerations' = lifelines_tables_1c_q_1_enumerations_loc)

table_and_columns[['2aq1']] <- list('location' = lifelines_tables_2a_q_1_results_loc,  'variables' =  variables_2aq1, 'enumerations' = lifelines_tables_2a_q_1_enumerations_loc)
table_and_columns[['2aq2']] <- list('location' = lifelines_tables_2a_q_2_results_loc,  'variables' =  variables_2aq2, 'enumerations' = lifelines_tables_2a_q_2_enumerations_loc)
table_and_columns[['2aqyouth']] <- list('location' = lifelines_tables_2a_q_youth_results_loc,  'variables' =  variables_2aqyouth, 'enumerations' = lifelines_tables_1a_q_youth_enumerations_loc)
table_and_columns[['2av1']] <- list('location' = lifelines_tables_2a_v_1_results_loc,  'variables' =  variables_2av1, 'enumerations' = lifelines_tables_2a_v_1_enumerations_loc)
table_and_columns[['2av2']] <- list('location' = lifelines_tables_2a_v_2_results_loc,  'variables' =  variables_2av2, 'enumerations' = lifelines_tables_2a_v_2_enumerations_loc)
table_and_columns[['2bq1']] <- list('location' = lifelines_tables_2b_q_1_results_loc,  'variables' =  variables_2bq1, 'enumerations' = lifelines_tables_2b_q_1_enumerations_loc)

table_and_columns[['3aq1']] <- list('location' = lifelines_tables_3a_q_1_results_loc,  'variables' =  variables_3aq1, 'enumerations' = lifelines_tables_3a_q_1_enumerations_loc)
table_and_columns[['3aq2']] <- list('location' = lifelines_tables_3a_q_2_results_loc,  'variables' =  variables_3aq2, 'enumerations' = lifelines_tables_3a_q_2_enumerations_loc)
table_and_columns[['3aqyouth']] <- list('location' = lifelines_tables_3a_q_youth_results_loc,  'variables' =  variables_3aqyouth, 'enumerations' = lifelines_tables_3a_q_youth_enumerations_loc)
# table_and_columns[['3av1']] <- list('location' = lifelines_tables_3a_v_1_results_loc,  'variables' =  variables_3av1, 'enumerations' = lifelines_tables_3a_v_1_enumerations_loc)
table_and_columns[['3bq1']] <- list('location' = lifelines_tables_3b_q_1_results_loc,  'variables' =  variables_3bq1, 'enumerations' = lifelines_tables_3b_q_1_enumerations_loc)

# get the metadata
metadata_lifelines <- get_relevant_metadata(table_and_columns, full_id_table_X1$psext)
# add the LL ID
metadata_lifelines <- cbind(data.frame(full_id_table_X1$ll), metadata_lifelines)

# get the table of STEMI mapping patients to controls
c_number_to_llid_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/stemi_anon_id_to_geno_id.tsv'
c_number_to_llid <- read.table(c_number_to_llid_loc, header = T, sep = ' ')

# subset to only the healthy controls
metadata_lifelines_controls <- metadata_lifelines[metadata_lifelines[['full_id_table_X1.ll']] %in% c_number_to_llid[['geno_id']], ]
# now also add the S-number
metadata_lifelines_controls <- cbind(c_number_to_llid[match(metadata_lifelines_controls[['full_id_table_X1.ll']], c_number_to_llid[['geno_id']]), c('sample'), drop = F], metadata_lifelines_controls)
# remove LL
metadata_lifelines_controls[['full_id_table_X1.ll']] <- NULL
# replace all the $6 with NA
for (col in colnames(metadata_lifelines_controls)) {
  metadata_lifelines_controls[[col]] <- gsub('\\$6', NA, metadata_lifelines_controls[[col]])
  metadata_lifelines_controls[[col]] <- gsub('\\$7', NA, metadata_lifelines_controls[[col]])
}
# make the neutrophil and lymphocyte info numeric
for (col in c('1av2.lymphoconc_result_all_m_1', '1av2.neutroconc_result_all_m_1', '2av2.lymphoconc_result_all_m_1', '2av2.neutroconc_result_all_m_1')) {
  metadata_lifelines_controls[[col]] <- as.numeric(metadata_lifelines_controls[[col]])
}
# write the full table
write.table(metadata_lifelines_controls, '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/stemi_control_lifeline_variables_wolf.tsv', sep = '\t', row.names = F, col.names = T)
