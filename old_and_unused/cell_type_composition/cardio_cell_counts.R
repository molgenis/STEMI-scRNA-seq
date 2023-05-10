library(ggplot2)
library(reshape2)
library(dplyr)

add_data <- function(data = NULL, table_path, condition, chemistry, celltypes_to_use = NULL) {
  cell_counts <- read.table(table_path, sep = "\t", header = T, row.names = 1)

  if(!is.null(celltypes_to_use)) cell_counts <- cell_counts[celltypes_to_use,]

  melted_counts <- melt(as.matrix(cell_counts), varnames = c("cell_type", "sample"), value.name = "count")
  melted_counts$condition <- condition
  melted_counts$chemistry <- chemistry

  if (is.null(data)) return(melted_counts)
  return(rbind(data, melted_counts))
}


plot_t_types <- function(cell_counts, CD4T_name='CD4T', CD8T_name='CD8T', plot_combined=T, log2transform=T){
  t_fraction_table <- NULL
  for(participant in unique(cell_counts$sample)){
    for(condition in unique(cell_counts[cell_counts$sample == participant, 'condition'])){
      # check if participant exists here
      if(nrow(cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, ]) > 0){
        CD4T_number <- 0
        CD8T_number <- 0
        # grab the values
        CD4T_number <- cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        CD8T_number <- cell_counts[cell_counts$cell_type == CD8T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        chem <- cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'chemistry']
        # calculate the T fraction
        t_fraction <- as.numeric(CD4T_number)/as.numeric(CD8T_number)
        # log2 transform if asked
        if(log2transform){
          t_fraction <- log2(t_fraction)
        }
        # add to table if exists, otherwise create the table
        if(is.null(t_fraction_table)){
          t_fraction_table <- data.frame(c(participant), c(t_fraction), c(condition), c(chem), stringsAsFactors = FALSE)
          colnames(t_fraction_table) <- c('sample', 'ratio', 'condition', 'chemistry')
        }
        else{
          t_fraction_table <- rbind(t_fraction_table, c(participant, t_fraction, condition, chem))
        }
      }
    }
  }
  # add the table as a whole once more, to get the combined V2+V3 ratios is requested
  if(plot_combined){
    # copy table
    t_fraction_table_both <- t_fraction_table
    # set chem to both
    t_fraction_table_both$chemistry <- 'both'
    # paste to original table
    t_fraction_table <- rbind(t_fraction_table, t_fraction_table_both)
  }
  # set y label based on log2 transformation
  ylabel <- 'ratio CD4T to CD8T'
  if(log2transform){
    ylabel <- 'log2 ratio CD4T to CD8T'
  }
  # do the actual plotting
  ggplot(t_fraction_table, aes(x=condition, y=as.numeric(ratio), fill=condition)) +
    geom_boxplot(color = "black", outlier.shape=NA,lwd=0.4, alpha=1) +
    #geom_boxplot(color = "black", lwd=0.4, alpha=1) +
    #scale_y_continuous(limits = c(0, 900)) +
    facet_grid(vars(chemistry)) +
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
          strip.text.y = element_text(colour = "black", size = 12, family = "Helvetica"),
          title = element_text(size = 12),
          axis.title.y = element_text(size = 12, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8, family = "Helvetica"),
          axis.text.x = element_blank(),
          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
    scale_fill_manual(values = c('blue', 'yellow', 'green'), name = "Timepoint") +
    labs(title = 'proportion of CD4T versus CD8T', y = ylabel)
  
}


plot_lymphoid_vs_myeloid <- function(cell_counts, plot_combined=T, log2transform=T, B_name='B', CD4T_name='CD4T', CD8T_name='CD8T', DC_name='DC', NK_name='NK', hemapoietic_stem_name='hemapoietic stem', megakaryocyte_name='megakaryocyte', monocyte_name='monocyte', plasma_B_name='plasma B', th17_name='Th17', tc17_name='Tc17'){
  fraction_table <- NULL
  for(participant in unique(cell_counts$sample)){
    for(condition in unique(cell_counts[cell_counts$sample == participant, 'condition'])){
      # check if participant exists here
      if(nrow(cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, ]) > 0){
        # grab the values
        B_number <- cell_counts[cell_counts$cell_type == B_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        CD4T_number <- cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        CD8T_number <- cell_counts[cell_counts$cell_type == CD8T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        DC_number <- cell_counts[cell_counts$cell_type == DC_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        NK_number <- cell_counts[cell_counts$cell_type == NK_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        hemapoietic_stemi_number <- cell_counts[cell_counts$cell_type == hemapoietic_stem_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        megakaryocyte_number <- cell_counts[cell_counts$cell_type == megakaryocyte_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        monocyte_number <- cell_counts[cell_counts$cell_type == monocyte_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        plasma_B_number <- cell_counts[cell_counts$cell_type == plasma_B_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        th17_number <- cell_counts[cell_counts$cell_type == th17_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        tc17_number <- cell_counts[cell_counts$cell_type == tc17_name & cell_counts$sample == participant & cell_counts$condition == condition, 'count']
        chem <- cell_counts[cell_counts$cell_type == CD4T_name & cell_counts$sample == participant & cell_counts$condition == condition, 'chemistry']
        # combine the myeloids
        myeloids <- c(monocyte_number, DC_number)
        # combine the lymphoids
        lymphoids <- c(B_number, plasma_B_number, CD4T_number, CD8T_number, NK_number, tc17_number, th17_number)
        # replace NAs with zero
        is.na(myeloids) <- 0
        is.na(lymphoids) <- 0
        # sum them
        myeloids_number <- sum(myeloids)
        lymphoids_number <- sum(lymphoids)
        fraction <- as.numeric(lymphoids_number)/as.numeric(myeloids_number)
        # log2 transform if asked
        if(log2transform){
          fraction <- log2(fraction)
        }
        # add to table if exists, otherwise create the table
        if(is.null(fraction_table)){
          fraction_table <- data.frame(c(participant), c(fraction), c(condition), c(chem), stringsAsFactors = FALSE)
          colnames(fraction_table) <- c('sample', 'ratio', 'condition', 'chemistry')
        }
        else{
          fraction_table <- rbind(fraction_table, c(participant, fraction, condition, chem))
        }
      }
    }
  }
  # add the table as a whole once more, to get the combined V2+V3 ratios is requested
  if(plot_combined){
    # copy table
    fraction_table_both <- fraction_table
    # set chem to both
    fraction_table_both$chemistry <- 'both'
    # paste to original table
    fraction_table <- rbind(fraction_table, fraction_table_both)
  }
  # set y label based on log2 transformation
  ylabel <- 'ratio lymphoid to myeloid'
  if(log2transform){
    ylabel <- 'log2 ratio lymphoid to myeloid'
  }
  # do the actual plotting
  ggplot(fraction_table, aes(x=condition, y=as.numeric(ratio), fill=condition)) +
    geom_boxplot(color = "black", outlier.shape=NA,lwd=0.4, alpha=1) +
    #geom_boxplot(color = "black", lwd=0.4, alpha=1) +
    #scale_y_continuous(limits = c(0, 900)) +
    facet_grid(vars(chemistry)) +
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
          strip.text.y = element_text(colour = "black", size = 12, family = "Helvetica"),
          title = element_text(size = 12),
          axis.title.y = element_text(size = 12, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8, family = "Helvetica"),
          axis.text.x = element_blank(),
          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
    scale_fill_manual(values = c('blue', 'yellow', 'green'), name = "Timepoint") +
    labs(title = 'proportion of lymphoid versus myeloid', y = ylabel)
        
}


metadata_to_cell_counts <- function(metadata, cell_type_column='cell_type', assignment.column='assignment.final', timepoint.column='timepoint.final'){
  formatted_counts <- NULL
  # check each condition
  for(condition in unique(metadata[[timepoint.column]])){
    # grab for this timepoint
    metadata_timepoint <- metadata[metadata[[timepoint.column]] == condition, ]
    # check each participant
    for(participant in unique(metadata_timepoint[[assignment.column]])){
      # grab for this participant
      metadata_tp_participant <- metadata_timepoint[metadata_timepoint[[assignment.column]] == participant, ]
      # check each cell type
      for(cell_type in unique(metadata_tp_participant[[cell_type_column]])){
        # finally grab the number we want
        cells_tp_part_ct <- nrow(metadata_tp_participant[metadata_tp_participant[[cell_type_column]] == cell_type, ])
        # only add if there are more than zero cells
        if(cells_tp_part_ct > 0){
          # grab the chemistry
          chem <- as.character(metadata_tp_participant$chem[1])
          # add to table
          if(is.null(formatted_counts)){
            formatted_counts <- data.frame(c(participant), c(cell_type), c(cells_tp_part_ct), c(condition), c(chem), stringsAsFactors = F)
            colnames(formatted_counts) <- c('sample', 'cell_type', 'count', 'condition', 'chemistry')
          }
          else{
            formatted_counts <- rbind(formatted_counts, c(participant, cell_type, cells_tp_part_ct, condition, chem))
          }
        }
      }
    }
  }
  return(formatted_counts)
}

metadata_to_cell_counts_colmethod <- function(metadata, cell_type_column='cell_type', assignment.column='assignment.final', timepoint.column='timepoint.final'){
  formatted_counts <- NULL
  # get all cell types
  possible_cell_types <- as.character(unique(metadata[[cell_type_column]]))
  # check each condition
  for(condition in unique(metadata[[timepoint.column]])){
    # grab for this timepoint
    metadata_timepoint <- metadata[metadata[[timepoint.column]] == condition, ]
    # check each participant
    for(participant in unique(metadata_timepoint[[assignment.column]])){
      # grab for this participant
      metadata_tp_participant <- metadata_timepoint[metadata_timepoint[[assignment.column]] == participant, ]
      # store the counts per cell type for the participant
      counts_part_and_condition <- c()
      # get the chemistry
      chem <- as.character(metadata_tp_participant$chem[1])
      # check each cell type
      for(cell_type in unique(metadata[[cell_type_column]])){
        # finally grab the number we want
        cells_tp_part_ct <- nrow(metadata_tp_participant[metadata_tp_participant[[cell_type_column]] == cell_type, ])
        # add this number of cells
        counts_part_and_condition <- c(counts_part_and_condition, cells_tp_part_ct)
      }
      # add the bulk number as well
      counts_part_and_condition <- c(counts_part_and_condition, nrow(metadata_tp_participant))
      # add to table
      if(is.null(formatted_counts)){
        formatted_counts <- data.frame(c(participant, counts_part_and_condition, condition, chem))
        # transpose
        formatted_counts <- t(formatted_counts)
        rownames(formatted_counts) <- NULL
        colnames(formatted_counts) <- c('sample', possible_cell_types,  'total', 'condition', 'chem')
      }
      else{
        formatted_counts <- rbind(formatted_counts, c(participant, counts_part_and_condition, condition, chem))
      }
    }
  }
  return(formatted_counts)
}



plot_mono1_vs_mono2 <- function(counts_data, cell_type_column='cell_type', mono1_name='mono 1', mono2_name='mono 2'){
  # subset to only have the mono 1 and mono 2 numbers
  counts_data_monos <- counts_data[counts_data[[cell_type_column]] == mono1_name | counts_data[[cell_type_column]] == mono2_name, ]
  ggplot(counts_data_monos, aes(x=cell_type, y=as.numeric(count), fill=cell_type)) +
    geom_boxplot(color = "black", outlier.shape=NA,lwd=0.4, alpha=1) +
    #geom_boxplot(color = "black", lwd=0.4, alpha=1) +
    scale_y_continuous(limits = c(0, 900)) +
    facet_grid(vars(chemistry), vars(condition)) +
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
          strip.text.y = element_text(colour = "black", size = 12, family = "Helvetica"),
          title = element_text(size = 12),
          axis.title.y = element_text(size = 12, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text.y = element_text(size = 8, family = "Helvetica"),
          axis.text.x = element_blank(),
          # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
    scale_fill_manual(values = c('brown', 'gray'), name = "Cell type")
  
}


features_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/'
features_loc <- '/data/cardiology/eQTL_mapping/features/'
#features_loc_v2 <- paste(features_loc, 'stemi_v2_full_res/', sep = '')
#features_loc_v3 <- paste(features_loc, 'stemi_v3_full_res/', sep = '')
features_loc_v2 <- paste(features_loc, 'stemi_v2_lowerres_20200625/', sep = '')
features_loc_v3 <- paste(features_loc, 'stemi_v3_lowerres_20200625/', sep = '')
conditions <- c('Baseline', 't24h', 't8w')


celltypes_to_use <- c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'Tc17', 'Th17', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B')

cell_counts <- NULL
for(condition in conditions){
  loc <- paste(features_loc_v2, condition, '/cell_counts.txt', sep = '')
  cell_counts <- add_data(data = cell_counts, table_path = loc, condition = condition, chemistry = "V2", celltypes_to_use = celltypes_to_use)
}
for(condition in conditions){
  loc <- paste(features_loc_v3, condition, '/cell_counts.txt', sep = '')
  cell_counts <- add_data(data = cell_counts, table_path = loc, condition = condition, chemistry = "V3", celltypes_to_use = celltypes_to_use)
}

cell_counts$condition = factor(cell_counts$condition, levels=conditions)
#colors = c("#153057", "#009ddb", "#e64b50", "#edba1b", "#71bc4b", "#965ec8", "#965ec8")
colors = c('red', 'blue', 'yellow', 'magenta', 'cyan', 'purple', 'green', 'orange', 'violet', 'tomato', 'pink')

##
## Boxplot, grid (condition x chemistry)
##
ggplot(cell_counts, aes(x=cell_type, y=count, fill=cell_type)) +
  geom_boxplot(color = "black", outlier.shape=NA,lwd=0.4, alpha=1) +
  #geom_boxplot(color = "black", lwd=0.4, alpha=1) +
  scale_y_continuous(limits = c(0, 900)) +
  facet_grid(vars(chemistry), vars(condition)) +
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
        strip.text.y = element_text(colour = "black", size = 12, family = "Helvetica"),
        title = element_text(size = 12),
        axis.title.y = element_text(size = 12, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "Helvetica"),
        axis.text.x = element_blank(),
       # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
  scale_fill_manual(values = colors, name = "Cell type")

ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/cell_counts_boxplot_condition_chemistry_20200622.png', dpi = 600, width = 10, height = 10)

##
## Boxplot, chemistry V2 vs V3
##
ggplot(cell_counts, aes(x=cell_type, y=count, fill=cell_type)) +
  geom_boxplot(color = "black", outlier.shape=NA,lwd=0.4, alpha=1) +
  #geom_boxplot(color = "black", lwd=0.4, alpha=1) +
  scale_y_continuous(limits = c(0, 900)) +
  facet_wrap(~chemistry) +
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 12, family = "Helvetica"),
        strip.text.y = element_text(colour = "black", size = 12, family = "Helvetica"),
        title = element_text(size = 12),
        axis.title.y = element_text(size = 12, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8, family = "Helvetica"),
        axis.text.x = element_blank(),
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "grey", fill=NA, size=1)) +
  scale_fill_manual(values = colors, name = "Cell type")

ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/cell_counts_chemistry_20200622.png', dpi = 600, width = 10, height = 10)



ggplot(filter(cell_counts, chemistry=="V3"), aes(sample, count)) +
  geom_col(aes(fill = cell_type), position="fill") +
  scale_fill_manual(values = colors, name = "Cell type") +
  facet_wrap(~condition) +
  ggtitle("V3 cell counts") +
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 12),
      strip.text.y = element_text(colour = "black", size = 12),
      title = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 8),
      panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5))

ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/cell_counts_v3_20200622.png', dpi = 600, width = 10, height = 10)



ggplot(filter(cell_counts, chemistry=="V2"), aes(sample, count)) +
  geom_col(aes(fill = cell_type), position="fill") +
  scale_fill_manual(values = colors, name = "Cell type") +
  facet_wrap(~condition) +
  ggtitle("V2 cell counts") +
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 12),
        strip.text.y = element_text(colour = "black", size = 12),
        title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, size = 6, vjust = 0.5))

ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/cell_counts_v2_20200622.png', dpi = 600, width = 10, height = 10)

ggplot(filter(cell_counts), aes(condition, count)) +
  geom_col(aes(fill = cell_type), position="fill") +
  scale_fill_manual(values = colors, name = "Cell type") +
  ggtitle("Cell counts V2 & V3") +
  facet_wrap(~chemistry, ncol = 1) +
  theme_minimal(base_family = "Helvetica") +
  theme(strip.text.x = element_text(colour = "black", size = 12),
        strip.text.y = element_text(colour = "black", size = 12),
        title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(size = 12))

ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/cell_counts_v2_and_v3_20200622.png', dpi = 600, width = 10, height = 10)


# grab the metadata
meta.data <- read.table('/data/cardiology/metadata/cardio.integrated_meta.data.tsv', sep='\t', header=T, row.names=1)
meta_data_counts <- metadata_to_cell_counts(meta.data)
plot_mono1_vs_mono2(meta_data_counts)
