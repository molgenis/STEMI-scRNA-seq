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

features_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/'
#features_loc_v2 <- paste(features_loc, 'stemi_v2_full_res/', sep = '')
#features_loc_v3 <- paste(features_loc, 'stemi_v3_full_res/', sep = '')
features_loc_v2 <- paste(features_loc, 'stemi_v2_lowerres_20200622/', sep = '')
features_loc_v3 <- paste(features_loc, 'stemi_v3_lowerres_20200622/', sep = '')
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



#cell_counts <- add_data(table_path = "../../1M_cells/data/cell_counts/v2/UT_cell_counts.txt", condition = "UT", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/3hCA_cell_counts.txt", condition = "3hCA", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/3hMTB_cell_counts.txt", condition = "3hMTB", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/3hPA_cell_counts.txt", condition = "3hPA", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/24hCA_cell_counts.txt", condition = "24hCA", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/24hMTB_cell_counts.txt", condition = "24hMTB", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v2/24hPA_cell_counts.txt", condition = "24hPA", chemistry = "V2", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/UT_cell_counts.txt", condition = "UT", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/3hCA_cell_counts.txt", condition = "3hCA", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/3hMTB_cell_counts.txt", condition = "3hMTB", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/3hPA_cell_counts.txt", condition = "3hPA", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/24hCA_cell_counts.txt", condition = "24hCA", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/24hMTB_cell_counts.txt", condition = "24hMTB", chemistry = "V3", celltypes_to_use = celltypes_to_use)
#cell_counts <- add_data(cell_counts, table_path = "../../1M_cells/data/cell_counts/v3/24hPA_cell_counts.txt", condition = "24hPA", chemistry = "V3", celltypes_to_use = celltypes_to_use)

# Order levels
#cell_counts$condition = factor(cell_counts$condition, levels=c("UT", "3hCA", "24hCA", "3hPA", "24hPA", "3hMTB", "24hMTB"))
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

# Investigate outlier
#cell_counts[cell_counts$count > 2000,]
#cell_counts[cell_counts$sample == "LLDeep_0471" & cell_counts$condition == "3hPA",]

#cell_counts_by_sample <- cell_counts %>% group_by(sample, cell_type) %>% summarize(total_count = sum(count))

#boxplot(filter(cell_counts_by_sample, cell_type=="monocyte")$total_count)

#ggplot(filter(cell_counts, condition=="UT" & chemistry=="V3"), aes(sample, count)) +
#  geom_col(aes(fill = cell_type), position="fill") +
#  scale_fill_manual(values = colors, name = "Cell type")

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
