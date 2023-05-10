######################
# libraries          #
######################

# none yet

####################
# Functions        #
####################

metadata_to_ggally_table <- function(metadata, cell_type_column='cell_type_lowerres', timepoint.column='timepoint.final', assignment.column='assignment.final', timepoints=c('UT', 'Baseline', 't24h', 't8w'), to_fraction=F){
  # get the unique cell types
  cell_types <- unique(metadata[[cell_type_column]])
  # get the unique participants
  participants <- unique(metadata[[assignment.column]])
  # get the unique timepoints
  timepoints <- unique(metadata[[timepoint.column]])
  # combinations of cell type and participant
  part_ct_comb <- paste(rep(participants, each = length(cell_types)), cell_types, sep = ".")
  # create matrix
  number_matrix <- matrix(, nrow = length(part_ct_comb), ncol=2+length(timepoints), dimnames = list(part_ct_comb, c('assignment', 'cell_type', timepoints)))
  # check the combinations
  for(cell_type in cell_types){
    for(participant in participants){
      # add to matrix
      number_matrix[paste(participant, cell_type, sep = '.'), 'assignment'] <- participant
      number_matrix[paste(participant, cell_type, sep = '.'), 'cell_type'] <- cell_type
      for(timepoint in timepoints){
        # get the number of cells
        nr_of_cells <- nrow(metadata[metadata[[cell_type_column]] == cell_type &
                                       metadata[[assignment.column]] == participant &
                                       metadata[[timepoint.column]] == timepoint,
                                     ])
        if(to_fraction){
          total_cells <- nrow(metadata[metadata[[assignment.column]] == participant &
                                         metadata[[timepoint.column]] == timepoint,
                                       ])
          nr_of_cells <- nr_of_cells/total_cells
        }
        # add to matrix
        number_matrix[paste(participant, cell_type, sep = '.'), timepoint] <- nr_of_cells
      }
    }
  }
  number_matrix <- data.frame(number_matrix)
  return(number_matrix)
}


add_peak_ck_mb <- function(clinvar_object){
  # subset to the ckmb variables
  ckmb <- clinvar[, grep('ck_mb', colnames(clinvar_object))]
  # change the NA variables to -1
  ckmb[is.na(ckmb)] <- -1
  # get the maximum
  max_ckmb <- apply(ckmb, 1, max)
  # replace a -1 with NA, this means there was no real variable
  max_ckmb[max_ckmb == -1] <- NA
  # paste onto the existing object
  clinvar_object$peak_ck_mb <- max_ckmb
  return(clinvar_object)
}

metadata_to_ct_table_per_column_per_part <- function(metadata, cell_type_column='cell_type_lowerres', timepoint.column='timepoint.final', assignment.column='assignment.final', timepoints_to_check=NULL, cell_types_to_check=NULL, na_to_zero=T, remove_zero_rows=T){
  # check either cell types provided or all
  cell_types <- cell_types_to_check
  if(is.null(cell_types)){
    cell_types <- unique(metadata[[cell_type_column]])
  }
  # check either timepoints provided or all
  timepoints <- timepoints_to_check
  if(is.null(timepoints)){
    timepoints <- unique(metadata[[timepoint.column]])
  }
  assignments <- unique(metadata[[assignment.column]])
  # create rownames
  rows <- paste(rep(timepoints, each = length(assignments)), assignments, sep = ".")
  # init the matrix
  numbers_table <- matrix(, nrow=length(rows), ncol=length(cell_types), dimnames = list(rows, cell_types))
  # check each timepoint
  for(timepoint in intersect(timepoints, unique(metadata[[timepoint.column]]))){
    # subset to that timepoint
    metadata.tp <- metadata[metadata[[timepoint.column]] == timepoint, ]
    # check each cell type
    for(cell_type in intersect(cell_types, unique(metadata.tp[[cell_type_column]]))){
      # subset to cell type
      metadata.tp.ct <- metadata.tp[metadata.tp[[cell_type_column]] == cell_type, ]
      # check each participant
      for(assignment in intersect(assignments, unique(metadata.tp.ct[[assignment.column]]))){
        # get the number of cells
        nr_of_cells <- nrow(metadata.tp.ct[metadata.tp.ct[[assignment.column]] == assignment, ])
        # create row name
        rowname <- paste(timepoint, assignment, sep = '.')
        # add to the matrix
        numbers_table[rowname, cell_type] <- nr_of_cells
      }
    }
  }
  # I like dataframes
  numbers_table <- data.frame(numbers_table)
  # na to zero if requested
  if(na_to_zero){
    numbers_table[is.na(numbers_table)] <- 0
  }
  # remove entries with no values if requested
  if(remove_zero_rows){
    # calculate sums
    sums_cells <- apply(numbers_table, 1, sum)
    # subset to non-zero
    numbers_table <- numbers_table[sums_cells > 0, ]
  }
  return(numbers_table)
}

per_column_ct_to_ggally <- function(ct_table_per_cond, sep.character='\\.'){
  # split the rows0
  conditions_and_participants <- strsplit(rownames(ct_table_per_cond) , sep.character)
  # set the vectors for conditions and participants
  conditions <- c()
  participants <- c()
  # grab the relevant items
  for(item in conditions_and_participants){
    # get the values
    cond <- item[1]
    part <- item[2]
    # add to vector
    conditions <- c(conditions, cond)
    participants <- c(participants, part)
  }
  # make unique
  conditions <- unique(conditions)
  participants <- unique(participants)
  # row names will be cell type and 
  rows <- paste(rep(participants, each = length(colnames(ct_table_per_cond))), colnames(ct_table_per_cond), sep = ".")
  # create table to store results
  numbers_table <- matrix(, nrow=length(rows), ncol=length(conditions)+2 , dimnames = list(rows, c('participant', 'cell_type', conditions)))
  # check each row
  for(row in rownames(ct_table_per_cond)){
    # get the condition and participant from the rowname
    cond_part <- strsplit(row, sep.character)
    cond <- cond_part[[1]][1]
    part <- cond_part[[1]][2]
    print(cond)
    print(part)
    # check each column, which is a cell type
    for(cell_type in colnames(ct_table_per_cond)){
      # grab the value
      ct_val <- ct_table_per_cond[row, cell_type]
      # add the values
      rowname <- paste(part, cell_type, sep='.')
      # set vars
      numbers_table[rowname, 'participant'] <- part
      numbers_table[rowname, 'cell_type'] <- cell_type
      numbers_table[rowname, cond] <- ct_val
    }
  }
  # to dataframe
  numbers_table <- data.frame(numbers_table)
  return(numbers_table)
}

####################
# Main Code        #
####################

# location of the metadata
meta.data.loc <- '/data/cardiology/metadata/cardio.integrated.20210301.metadata.tsv'
# location of the clinical variables
clinvar.loc <- '/data/cardiology/metadata/SingleCellSequencing_20210114.tsv'

# read into table
meta.data <- read.table(meta.data.loc, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
# use the gally method to plot the cell numbers
cell_numbers_gally_stemi <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ])
# do the same, but with fractions
cell_numbers_gally_stemi_frac <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ], to_fraction = T)
# now with the lfc of fractions
cell_numbers_gally_stemi_frac_baselinescaled <- cell_numbers_gally_stemi_frac[, c('assignment', 'cell_type')]
cell_numbers_gally_stemi_frac_baselinescaled$t24h_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['Baseline']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled$t24h_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['t24h']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled$t8w_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['t8w']) / as.numeric(x['Baseline']))})
# high resolution as well
cell_numbers_gally_stemi_frac_hr <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ], to_fraction = T, cell_type_column = 'cell_type')
cell_numbers_gally_stemi_frac_baselinescaled_hr <- cell_numbers_gally_stemi_frac_hr[, c('assignment', 'cell_type')]
cell_numbers_gally_stemi_frac_baselinescaled_hr$t24h_t0 <- apply(cell_numbers_gally_stemi_frac_hr, 1, function(x){log2(as.numeric(x['Baseline']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled_hr$t24h_t0 <- apply(cell_numbers_gally_stemi_frac_hr, 1, function(x){log2(as.numeric(x['t24h']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled_hr$t8w_t0 <- apply(cell_numbers_gally_stemi_frac_hr, 1, function(x){log2(as.numeric(x['t8w']) / as.numeric(x['Baseline']))})


# read the clinvar table
clinvar <- read.table(clinvar.loc, header = T, sep = '\t', dec = '.')
# calculate the peak ck_mb
clinvar <- add_peak_ck_mb(clinvar)
# add log of peak and baseline ckmb
clinvar$log_peak_ck_mb <- log(clinvar$peak_ck_mb)
clinvar$log_ck_mb <- log(clinvar$ck_mb)
# calculate the age
clinvar$date_birth <- ymd(clinvar$date_birth)
clinvar$admission_date <- ymd(clinvar$admission_date)
clinvar$age <- interval(start= clinvar$date_birth, end=clinvar$admission_date)/duration(n=1, unit="years")
# merge the cell numbers with the clinical variables
clinvar_wfracs <- merge(clinvar, cell_numbers_gally_stemi_frac_baselinescaled, by.x = 'record_id', by.y = 'assignment', all = T)
summary(lm(formula = t24h_t0 ~ log_peak_ck_mb+age+gender+ischemia_time+log_ck_mb, data=clinvar_wfracs[!is.na(clinvar_wfracs$t24h_t0) & clinvar_wfracs$cell_type == 'monocyte', ]))
summary(lm(formula = t24h_t0 ~ log_peak_ck_mb+age+gender+log_ck_mb, data=clinvar_wfracs[!is.na(clinvar_wfracs$t24h_t0) & clinvar_wfracs$cell_type == 'monocyte', ]))
summary(lm(formula = t24h_t0 ~ peak_ck_mb+age+gender+ck_mb, data=clinvar_wfracs[!is.na(clinvar_wfracs$t24h_t0) & clinvar_wfracs$cell_type == 'monocyte', ]))
# now to the same stuff with hr and cmono
clinvar_wfracs_hr <- merge(clinvar, cell_numbers_gally_stemi_frac_baselinescaled_hr, by.x = 'record_id', by.y = 'assignment', all = T)
summary(lm(formula = t24h_t0 ~ log_peak_ck_mb+age+gender+ischemia_time+log_ck_mb, data=clinvar_wfracs_hr[!is.na(clinvar_wfracs_hr$t24h_t0) & clinvar_wfracs_hr$cell_type == 'cMono', ]))

