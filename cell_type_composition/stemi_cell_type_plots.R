####################
# libraries        #
####################

library(ggplot2)
library(reshape2)
library(dplyr)

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
      total_cells <- nrow(metadata.tp)
      # check each cell type
      for(cell_type in unique(metadata.tp.part[[cell_type_column]])){
        # get that number of cells
        nr_cells <- nrow(metadata.tp.part[metadata.tp.part[[cell_type_column]] == cell_type, ])
        # turn into that dataframe
        nr_cells_df <- data.frame(participant=c(participant), condition=c(timepoint), cell_type=c(cell_type), number=c(nr_cells), fraction=c(nr_cells/total_cells))
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
          scaled_df <- data.frame(participant=c(participant), condition=c(timepoint_to_scale), cell_type=c(cell_type), number=c(scaled))
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

plot_ct_numbers_boxplot <- function(numbers_table, conditions_to_plot=c('UT', 'Baseline', 't24h', 't8w'), cell_types_to_plot=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, to_fraction=F, pointless=F, legendless=F){
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
  # fetch colors
  cc <- get_color_coding_dict()
  # set colors based on condition
  colScale <- scale_fill_manual(name = "condition",values = unlist(cc[numbers_table$condition]))
  # create the plot
  p <- ggplot(data=numbers_table, aes(x=condition, y=number, fill=condition)) + 
    geom_boxplot() + 
    colScale +
    geom_jitter(size = 0.5, alpha = 0.5) +
    facet_grid(. ~ cell_type) +
    ylab(ylabel)
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
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
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
  label_dict[['UT']] <- 'HC'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
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
  return(label_dict)
}

####################
# Main Code        #
####################

# location of the metadata
meta.data.loc <- '/data/cardiology/metadata/cardio.integrated.20210301.metadata.tsv'
# read into table
meta.data <- read.table(meta.data.loc, sep = '\t', header = T, row.names = 1)
# get the cell numbers
cell_numbers <- metadata_to_ct_table(meta.data)
# plot the cell numbers
plot_ct_numbers_boxplot(numbers_table = cell_numbers, to_fraction = F, legendless = F, pointless = T)
# use the gally method to plot the cell numbers
cell_numbers_gally_stemi <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ])
ggparcoord(cell_numbers_gally_stemi[cell_numbers_gally_stemi$cell_type == 'monocyte', ], columns = c(4,3,5), groupColumn = 'assignment')
# do the same, but with fractions
cell_numbers_gally_stemi_frac <- metadata_to_ggally_table(meta.data[meta.data$orig.ident == 'stemi_v2' | meta.data$orig.ident == 'stemi_v3', ], to_fraction = T)
ggparcoord(cell_numbers_gally_stemi_frac[cell_numbers_gally_stemi_frac$cell_type == 'monocyte', ], columns = c(4,3,5), groupColumn = 'assignment')
# change to scaled cell numbers
cell_numbers_baselinescaled <- scale_cell_numbers_to_condition(cell_numbers, c('t24h', 't8w'), 'Baseline')
# plot these scaled numbers
plot_ct_numbers_boxplot(numbers_table = cell_numbers_baselinescaled, pointless = T)


