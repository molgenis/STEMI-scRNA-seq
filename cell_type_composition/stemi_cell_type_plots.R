####################
# libraries        #
####################

library(ggplot2)
library(GGally)

####################
# Functions        #
####################

metadata_to_ct_table <- function(metadata, cell_type_column='cell_type_lowerres', timepoint.column='timepoint.final', assignment.column='assignment.final', timepoints=c('UT', 'Baseline', 't24h', 't8w')){
  # init table
  cell_type_table <- NULL
  # check each timepoint
  for(timepoint in intersect(timepoints, unique(metadata[[timepoint.column]]))){
    # subset to that timepoint
    metadata.tp <- metadata[metadata[[timepoint.column]] == timepoint, ]
    # check each participant
    for(participant in unique(metadata.tp[[assignment.column]])){
      # subset to that participant
      metadata.tp.part <- metadata.tp[metadata.tp[[assignment.column]] == participant, ]
      # get the total number of cells
      total_cells <- nrow(metadata.tp.part)
      # check each cell type
      for(cell_type in unique(metadata.tp.part[[cell_type_column]])){
        # get that number of cells
        nr_cells <- nrow(metadata.tp.part[metadata.tp.part[[cell_type_column]] == cell_type, ])
        # turn into that dataframe
        nr_cells_df <- data.frame(participant=c(participant), condition=c(timepoint), cell_type=c(cell_type), number=c(nr_cells), fraction=c(nr_cells/total_cells), stringsAsFactors = F)
        # add to existing dataframe
        if(is.null(cell_type_table)){
          cell_type_table <- nr_cells_df
        }
        else{
          cell_type_table <- rbind(cell_type_table, nr_cells_df)
        }
      }
    }
  }
  return(cell_type_table)
}

scale_cell_numbers_to_condition <- function(cell_numbers, timepoints_to_scale, timepoint_to_scale_by){
  # init table
  scaled_number_table <- NULL
  # check each cell type
  for(cell_type in unique(cell_numbers$cell_type)){
    # subset to cell type
    cell_numbers_ct <- cell_numbers[cell_numbers$cell_type == cell_type, ]
    # check each participant
    for(participant in unique(cell_numbers_ct$participant)){
      # subset to participant
      cell_numbers_ct_part <- cell_numbers_ct[cell_numbers_ct$participant == participant, ]
      # check each timepoint to scale
      for(timepoint_to_scale in timepoints_to_scale){
        # check we have what to scale and what to scale by
        if(nrow(cell_numbers_ct_part[cell_numbers_ct_part$condition == timepoint_to_scale_by, ]) > 0 & nrow(cell_numbers_ct_part[cell_numbers_ct_part$condition == timepoint_to_scale, ]) > 0){
          # get the scaled value
          scaled <- cell_numbers_ct_part[cell_numbers_ct_part$condition == timepoint_to_scale, 'number'] / cell_numbers_ct_part[cell_numbers_ct_part$condition == timepoint_to_scale_by, 'number']
          # log2 transform
          scaled <- log2(scaled)
          # add to table
          scaled_df <- data.frame(participant=c(participant), condition=c(timepoint_to_scale), cell_type=c(cell_type), number=c(scaled), stringsAsFactors = F)
          if(is.null(scaled_number_table)){
            scaled_number_table <- scaled_df
          }
          else{
            scaled_number_table <- rbind(scaled_number_table, scaled_df)
          }
        }
      }
    }
  }
  return(scaled_number_table)
}

plot_ct_numbers_boxplot <- function(numbers_table, conditions_to_plot=c('UT', 'Baseline', 't24h', 't8w'), cell_types_to_plot=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, to_fraction=F, pointless=F, legendless=F, to_pct=F, ylim=NULL){
  # subset to want we want to plot
  numbers_table <- numbers_table[numbers_table$condition %in% conditions_to_plot, ]
  numbers_table <- numbers_table[numbers_table$cell_type %in% cell_types_to_plot, ]
  # use prettier labels if requested
  if(use_label_dict){
    numbers_table$condition <- as.vector(unlist(label_dict()[numbers_table$condition]))
    numbers_table$cell_type <- as.vector(unlist(label_dict()[numbers_table$cell_type]))
  }
  # set y label
  ylabel <- 'number of cells'
  # use fractions if requested
  if(to_fraction){
    # fraction is already in there
    numbers_table$number <- numbers_table$fraction
    # ovewrite label as well
    ylabel <- 'fraction of cells'
  }
  if(to_pct){
    # based on fraction
    numbers_table$number <- numbers_table$fraction*100
    # ovewrite label as well
    ylabel <- 'percentage of cells'
  }
  # fetch colors
  cc <- get_color_coding_dict()
  # set colors based on condition
  colScale <- scale_fill_manual(name = "condition",values = unlist(cc[numbers_table$condition]))
  # create the plot
  p <- ggplot(data=numbers_table, aes(x=condition, y=number, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + 
    colScale +
    geom_jitter(size = 0.5, alpha = 0.5) +
    facet_grid(. ~ cell_type) +
    ylab(ylabel)
  # add xlimit if requested
  if(!is.null(ylim)){
    p <- p + ylim(ylim)
  }
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

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

get_cell_types_per_condition <- function(numbers_table){
  aggregate_df <- NULL
  # check per condition
  for(condition in unique(numbers_table$condition)){
    # subset to condition
    numbers_condition <- numbers_table[numbers_table$condition == condition, ]
    # count total number of cells in condition
    nr_condition <- sum(numbers_condition$number)
    # check cell types
    for(cell_type in unique(numbers_condition$cell_type)){
      # subset to cell type
      numbers_condition_cell_type <- numbers_condition[numbers_condition$cell_type == cell_type, ]
      # count total number of cells of cell type in condition
      nr_condition_cell_type <- sum(numbers_condition_cell_type$number)
      # calculate the fraction
      fraction <- nr_condition_cell_type / nr_condition
      # add to dataframe
      aggregate_row <- data.frame(condition=c(condition), cell_type=c(cell_type), number=c(nr_condition_cell_type), fraction=c(fraction))
      if(is.null(aggregate_df)){
        aggregate_df <- aggregate_row
      }
      else{
        aggregate_df <- rbind(aggregate_df, aggregate_row)
      }
    }
  }
  return(aggregate_df)
}

plot_cell_type_per_condition_bars <- function(numbers_per_cond, to_fraction=T, use_label_dict=T, split_chem=F){
  # get prettier labels if requested
  if(use_label_dict){
    numbers_per_cond$condition <- as.vector(unlist(label_dict()[numbers_per_cond$condition]))
    numbers_per_cond$cell_type <- as.vector(unlist(label_dict()[numbers_per_cond$cell_type]))
  }
  # init plot
  p <- NULL
  # use either the absolute or relative numbers
  if(to_fraction){
    p <- ggplot(numbers_per_cond, aes(x=condition, y=fraction, fill=cell_type))
  }
  else{
    p <- ggplot(numbers_per_cond, aes(x=condition, y=number, fill=cell_type))
  }
  # set colors based on condition
  cc <- get_color_coding_dict()
  colScale <- scale_fill_manual(name = 'cell type',values = unlist(cc[numbers_per_cond$cell_type]))
  # create plot
  p <- p + geom_bar(position='stack', stat='identity') +
    colScale
  # split by chem if possible
  if(split_chem){
    p <- p + facet_grid(. ~ chem)
  }
  return(p)
}

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

metadata_to_ct_table_per_column <- function(metadata, cell_type_column='cell_type_lowerres', timepoint.column='timepoint.final', timepoints_to_check=NULL, cell_types_to_check=NULL){
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
  # init the matrix
  numbers_table <- matrix(, nrow=length(timepoints), ncol=length(cell_types), dimnames = list(timepoints, cell_types))
  # check each timepoint
  for(timepoint in intersect(timepoints, unique(metadata[[timepoint.column]]))){
    # subset to that timepoint
    metadata.tp <- metadata[metadata[[timepoint.column]] == timepoint, ]
    # check each cell type
    for(cell_type in intersect(cell_types, unique(metadata.tp[[cell_type_column]]))){
      # get the number of cells
      nr_of_cells <- nrow(metadata.tp[metadata.tp[[cell_type_column]] == cell_type, ])
      # add to the matrix
      numbers_table[timepoint, cell_type] <- nr_of_cells
    }
  }
  # I like dataframes
  numbers_table <- data.frame(numbers_table)
  return(numbers_table)
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


reduce_table <- function(table_from, columns_to_keep=c('timepoint', 'condition', 'cell_type'), reducing_column='number'){
  # get the column index we are looking at
  column_index_to_keep <- length(columns_to_keep)-1
  # always grab the last column
  relevant_column = columns_to_keep[column_index_to_keep]
  # init table
  aggregate_table <- NULL
  # check each possible variable
  for(variable in unique(table_from[[relevant_column]])){
    # get the relevant columns for this one
    subset_table <- table_from[table_from[[relevant_column]] == variable, , drop=F]
    
    aggregate_rows <- NULL
    # if we won't nest deeper, we can actually sum this
    if(column_index_to_keep == 0){
      # calculate the sum
      aggregated_number <- sum(subset_table[[reducing_column]])
      # make a row
      aggregate_rows <- data.frame(number=c(aggregated_number))
    }
    else{
      # if there are still columns to 
      aggregated_rows <- reduce_table(subset_table, columns_to_keep[0:column_index_to_keep], reducing_column)
    }
    # add to aggregate table
    if(is.null(aggregate_table)){
      aggregate_table <- aggregate_rows
    }
    else{
      aggregate_table <- rbind(aggregate_table, aggregate_rows)
    }
    aggregate_table[[relevant_column]] <- variable
  }
  return(aggregate_table)
}

get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "paleturquoise3"
  color_coding[["t0-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

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
  label_dict[['UT']] <- 'HC'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
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

####################
# Main Code        #
####################

# location of the metadata
meta.data.loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/cardio.integrated.20210301.metadata.tsv'
# read into table
meta.data <- read.table(meta.data.loc, sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
# get the cell numbers
cell_numbers <- metadata_to_ct_table(meta.data)
# plot the cell numbers
plot_ct_numbers_boxplot(numbers_table = cell_numbers, legendless = F, pointless = T, to_pct = T, ylim=c(0,100))
# change to scaled cell numbers
cell_numbers_baselinescaled <- scale_cell_numbers_to_condition(cell_numbers, c('t24h', 't8w'), 'Baseline')
# plot these scaled numbers
plot_ct_numbers_boxplot(numbers_table = cell_numbers_baselinescaled, pointless = T)
# for higher res as well
cell_numbers_hr <- metadata_to_ct_table(meta.data, cell_type_column = 'cell_type')
plot_ct_numbers_boxplot(numbers_table = cell_numbers_hr[cell_numbers_hr$cell_type %in% c('cMono', 'NKdim'), ], to_fraction = T, legendless = F, pointless = T)


# use the gally method to plot the cell numbers
cell_numbers_gally_stemi <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ])
ggparcoord(cell_numbers_gally_stemi[cell_numbers_gally_stemi$cell_type == 'monocyte', ], columns = c(4,3,5), groupColumn = 'assignment')
# do the same, but with fractions
cell_numbers_gally_stemi_frac <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ], to_fraction = T)
ggparcoord(cell_numbers_gally_stemi_frac[cell_numbers_gally_stemi_frac$cell_type == 'monocyte', ], columns = c(4,3,5), groupColumn = 'assignment')
# now with the lfc of fractions
cell_numbers_gally_stemi_frac_baselinescaled <- cell_numbers_gally_stemi_frac[, c('assignment', 'cell_type')]
#cell_numbers_gally_stemi_frac_baselinescaled$t0_t0 <- log2(as.numeric(cell_numbers_gally_stemi_frac$Baseline) / as.numeric(cell_numbers_gally_stemi_frac$Baseline))
#cell_numbers_gally_stemi_frac_baselinescaled$t24h_t0 <- log2(as.numeric(cell_numbers_gally_stemi_frac$t24h) / as.numeric(cell_numbers_gally_stemi_frac$Baseline))
#cell_numbers_gally_stemi_frac_baselinescaled$t8w_t0 <- log2(as.numeric(cell_numbers_gally_stemi_frac$t8w) / as.numeric(cell_numbers_gally_stemi_frac$Baseline))
cell_numbers_gally_stemi_frac_baselinescaled$t24h_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['Baseline']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled$t24h_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['t24h']) / as.numeric(x['Baseline']))})
cell_numbers_gally_stemi_frac_baselinescaled$t8w_t0 <- apply(cell_numbers_gally_stemi_frac, 1, function(x){log2(as.numeric(x['t8w']) / as.numeric(x['Baseline']))})
ggparcoord(cell_numbers_gally_stemi_frac_baselinescaled[cell_numbers_gally_stemi_frac_baselinescaled$cell_type == 'monocyte' & cell_numbers_gally_stemi_frac_baselinescaled$assignment != 'TEST_81', ], columns = c(3,4,5), groupColumn = 'assignment', scale = 'globalminmax')


# set the fraction as the number
cell_numbers_fracasnr <- cell_numbers
cell_numbers_fracasnr$number <- cell_numbers_fracasnr$fraction
# scale that
cell_numbers_fracasnr_baselinescaled <- scale_cell_numbers_to_condition(cell_numbers_fracasnr, c('t24h', 't8w'), 'Baseline')
# set label correctly again
cell_numbers_fracasnr_baselinescaled$fraction <- cell_numbers_fracasnr_baselinescaled$number
# plot
plot_ct_numbers_boxplot(numbers_table = cell_numbers_fracasnr_baselinescaled, pointless = T, to_fraction = T)

# get cell numbers regardless of the participant
numbers_per_cond <- get_cell_types_per_condition(numbers_table = cell_numbers)
plot_cell_type_per_condition_bars(numbers_per_cond)
# now for v2 and v3 separately
cell_numbers_v2 <- metadata_to_ct_table(meta.data[meta.data$chem == 'V2', ])
cell_numbers_v3 <- metadata_to_ct_table(meta.data[meta.data$chem == 'V3', ])
numbers_per_cond_v2 <- get_cell_types_per_condition(numbers_table = cell_numbers_v2)
numbers_per_cond_v3 <- get_cell_types_per_condition(numbers_table = cell_numbers_v3)
numbers_per_cond_v2$chem <- 'V2'
numbers_per_cond_v3$chem <- 'V3'
numbers_per_cond_chemsplit <- rbind(numbers_per_cond_v2, numbers_per_cond_v3)
plot_cell_type_per_condition_bars(numbers_per_cond_chemsplit, split_chem = T)
ggsave('~/Desktop/cardio.integrated.20210301.cell_type_proportions_lowerres_overall_per_chem.pdf', width = 10, height = 5)

# create the cell type number tables as required in the supplements
ct_tbl <- metadata_to_ct_table_per_column(meta.data)
ct_tbl_hr <- metadata_to_ct_table_per_column(meta.data, cell_type_column = 'cell_type')
# make the col and rownames prettier
colnames(ct_tbl) <- as.vector(unlist(label_dict()[colnames(ct_tbl)]))
rownames(ct_tbl) <- as.vector(unlist(label_dict()[rownames(ct_tbl)]))
colnames(ct_tbl_hr) <- as.vector(unlist(label_dict()[colnames(ct_tbl_hr)]))
rownames(ct_tbl_hr) <- as.vector(unlist(label_dict()[rownames(ct_tbl_hr)]))

model_ncMono_t24ht0 <- lm(data=cell_numbers_gally_stemi_frac_baselinescaled_hr_ncmono_clinvar[!is.na(cell_numbers_gally_stemi_frac_baselinescaled_hr_ncmono_clinvar$t24h_t0), ], formula=t24h_t0 ~ logpeakckmb+age+gender+ck_mb)

# get the unique lane and participant combinations
participant_condition_lane <- unique(meta.data[, c('lane', 'assignment.final', 'timepoint.final')])
# add the lane
cell_numbers[['lane']] <- apply(cell_numbers, 1, FUN = function(x){
  lane <- participant_condition_lane[participant_condition_lane[['assignment.final']] == x['participant'] & participant_condition_lane[['timepoint.final']] == x['condition'], 'lane'][1]
  return(lane)
})
cell_numbers_hr[['lane']] <- apply(cell_numbers_hr, 1, FUN = function(x){
  lane <- participant_condition_lane[participant_condition_lane[['assignment.final']] == x['participant'] & participant_condition_lane[['timepoint.final']] == x['condition'], 'lane'][1]
  return(lane)
})
# add major/minor label
cell_numbers[['level']] <- 'major'
cell_numbers_hr[['level']] <- 'minor'
# turn into one table
cell_numbers_both <- rbind(cell_numbers, cell_numbers_hr)
cell_numbers_both[['condition']] <- as.vector(unlist(list('UT' = 'control', 'Baseline' = 't0', 't24h' = 't24h', 't8w' = 't6-8w')[cell_numbers_both[['condition']]]))
# write result
write.table(cell_numbers_both[, c('level', 'lane', 'condition', 'participant', 'cell_type', 'number')], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/stemi_cell_numbers_overview.tsv', sep = '\t', row.names = F, col.names = T)
