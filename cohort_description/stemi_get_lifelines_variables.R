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
      print(table)
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
ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_pseudo_ids.txt'
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/rsc01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/rsc01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
pseudo_ext_to_cytoid_map_loc <- '/groups/umcg-lifelines/rsc01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat'

# location of the questionaire data
lifelines_tables_loc <- '/groups/umcg-lifelines/rsc01/releases/pheno_lifelines/v2/results/'
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
lifelines_tables_loc <- '/groups/umcg-lifelines/rsc01/releases/pheno_lifelines/v2/enumerations/'
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


# read the tables
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_int_to_ugli_map <- read.table(pseudo_int_to_ugli_map_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_ext_to_cytoid_map <- read.table(pseudo_ext_to_cytoid_map_loc, sep = '\t', header = T, stringsAsFactors = F)

# combine them
full_id_table <- ll_to_pseudo_ext(ll_to_pseudo, pseudo_int_to_ext)

# set the appended LLDeep names
full_id_table_X1 <- full_id_table
full_id_table_X1$ll <- paste('X1_', full_id_table_X1$ll, sep = '')

# add the UGLI identifiers as well
full_id_table_X1 <- pseudo_int_to_ugli(full_id_table_X1, pseudo_int_to_ugli_map)

# add the cyto ID identifiers as well
full_id_table_X1 <- pseudo_ext_to_cytoid(full_id_table_X1, pseudo_ext_to_cytoid_map)

# we want to add the exp and age as well
exp_to_ll_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
age_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_age_gender_expidonly.tsv'

# read the files
exp_to_ll <- read.table(exp_to_ll_loc, sep = '\t', header = T, stringsAsFactors = F)
exp_to_ll$ll <- paste('X1_', exp_to_ll$LLD.ID, sep = '')
age_metadata <- read.table(age_metadata_loc, sep = '\t', header = T, stringsAsFactors = F)

# link the exp id to the data
full_id_table_X1$exp <- exp_to_ll[match(full_id_table_X1$ll, exp_to_ll$ll), 'ExpNr']
# link the age and gender to this
full_id_table_X1 <- cbind(full_id_table_X1, age_metadata[match(full_id_table_X1$exp, age_metadata$ExpNr), c('Gender', 'age_range')])

# check where we have lanes that are complete
stim_mapping_loc <- "/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_lane_to_tp.tsv"
# check_missingnes(stim_mapping_loc, full_id_table_X1[!is.na(full_id_table_X1$ugli), ]$exp)


####################
# get the ll data  #
####################

# the table that has the names of the variable I want
variable_mapping_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/stemi_lifelines_variables_codes.tsv'
variable_mapping <- read.table(variable_mapping_loc, sep = '\t', header = F)
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
stemi_to_control_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/stemi_age_gender_match_wtestid.tsv'
stemi_to_control <- read.table(stemi_to_control_loc, header = T, sep = '\t')

# subset to only the healthy controls
metadata_lifelines_controls <- metadata_lifelines[metadata_lifelines[['full_id_table_X1.ll']] %in% paste('X1_', stemi_to_control[['ll_match_id']], sep = ''), ]
# calculate BMI
metadata_lifelines_controls[['BMI']] <- as.numeric(metadata_lifelines_controls[['2av1.bodyweight_kg_all_m_1']]) / ((as.numeric(metadata_lifelines_controls[['2av1.bodylength_cm_all_m_1']])/100)^2)
mean(metadata_lifelines_controls[['BMI']][!is.na(metadata_lifelines_controls[['BMI']])])
sd(metadata_lifelines_controls[['BMI']][!is.na(metadata_lifelines_controls[['BMI']])])
# write the full table
write.table(metadata_lifelines_controls, '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/stemi_control_lifeline_variables.tsv', sep = '\t', row.names = F, col.names = T)

# calculate hypertension
sum(!is.na(metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']]) & metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']] == 'yes')
sum(!is.na(metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']]) & metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']] == 'yes') / sum(!is.na(metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']]) & (metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']] == 'yes' | metadata_lifelines_controls[['3aq1.hypertension_presence_adu_q_1']] == 'no'))

# calculate Hypercholesterolemia
sum(!is.na(metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']]) & metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes')
sum(!is.na(metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']]) & metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes') / sum(!is.na(metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']]) & (metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes' | metadata_lifelines_controls[['1aq1.highcholesterol_presence_adu_q_1']] == 'no'))

# check how many are smoking
sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) >= as.numeric(metadata_lifelines_controls[['1aq2.age']]))
# check how many stopped smoking
sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) < as.numeric(metadata_lifelines_controls[['1aq2.age']]))
# check how many never smoked
sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']] == '$6')

# check how many stopped smoking
sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) < as.numeric(metadata_lifelines_controls[['1aq2.age']]))
sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) < as.numeric(metadata_lifelines_controls[['1aq2.age']])) /
  sum(!is.na(metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_startage_adu_q_1_01']] != '$6' &
        !is.na(metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
        !is.na(metadata_lifelines_controls[['1aq2.age']]))

# get the metadata
cardio_integrated_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/cardio.integrated.20210301.metadata.tsv'
cardio_integrated_metadata <- read.table(cardio_integrated_metadata_loc, header = T, sep = '\t', stringsAsFactors = F)
# get the chemistry per participant
part_chem <- unique(cardio_integrated_metadata[, c('assignment.final', 'chem')])
# clean the table up
rownames(part_chem) <- NULL
colnames(part_chem) <- c('participant', 'chemistry')
# get the number of females in the V2 controls
full_id_table_X1_v2 <- full_id_table_X1[full_id_table_X1$ll %in% paste('X1_', part_chem[part_chem[['chemistry']] == 'V2', 'participant'], sep = ''), ]
# get number of females
nrow(full_id_table_X1_v2[full_id_table_X1_v2[['Gender']] == 'F', ])
nrow(full_id_table_X1_v2[full_id_table_X1_v2[['Gender']] == 'F', ]) / nrow(full_id_table_X1_v2)
# and age
exp_to_ll_v2 <- exp_to_ll[exp_to_ll$LLD.ID %in% part_chem[part_chem[['chemistry']] == 'V2', 'participant'], ]
mean(exp_to_ll_v2[['Age']])
sd(exp_to_ll_v2[['Age']])

# subset to v2
metadata_lifelines_controls_v2 <- metadata_lifelines_controls[metadata_lifelines_controls$full_id_table_X1.ll %in% paste('X1_', part_chem[part_chem[['chemistry']] == 'V2', 'participant'], sep = ''), ]
# calculate BMI
metadata_lifelines_controls_v2[['BMI']] <- as.numeric(metadata_lifelines_controls_v2[['2av1.bodyweight_kg_all_m_1']]) / ((as.numeric(metadata_lifelines_controls_v2[['2av1.bodylength_cm_all_m_1']])/100)^2)
mean(metadata_lifelines_controls_v2[['BMI']][!is.na(metadata_lifelines_controls_v2[['BMI']])])
sd(metadata_lifelines_controls_v2[['BMI']][!is.na(metadata_lifelines_controls_v2[['BMI']])])
# calculate hypertension
sum(!is.na(metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']]) & metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']] == 'yes')
sum(!is.na(metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']]) & metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']] == 'yes') / sum(!is.na(metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']]) & (metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']] == 'yes' | metadata_lifelines_controls_v2[['3aq1.hypertension_presence_adu_q_1']] == 'no'))
# calculate Hypercholesterolemia
sum(!is.na(metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']]) & metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes')
sum(!is.na(metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']]) & metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes') / sum(!is.na(metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']]) & (metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']] == 'yes' | metadata_lifelines_controls_v2[['1aq1.highcholesterol_presence_adu_q_1']] == 'no'))
# current smoker (smoking startage is not missing and the end age is not smaller than metadata_lifelines_controls_v2 current age)
sum(!is.na(metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls_v2[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']]) >= as.numeric(metadata_lifelines_controls_v2[['1aq2.age']]))
# check how many stopped smoking
sum(!is.na(metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']]) & metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']] != '$6' &
      !is.na(metadata_lifelines_controls_v2[['1aq2.age']]) &
      as.numeric(metadata_lifelines_controls_v2[['1aq2.smoking_endage_adu_qc_1_01']]) < as.numeric(metadata_lifelines_controls_v2[['1aq2.age']]))
# check how many never smoked
sum(!is.na(metadata_lifelines_controls_v2[['1aq2.smoking_startage_adu_q_1_01']]) & metadata_lifelines_controls_v2[['1aq2.smoking_startage_adu_q_1_01']] == '$6')

