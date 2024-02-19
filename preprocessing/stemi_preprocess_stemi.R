#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_preprocess_stemi.R
# Function: preprocess STEMI samples into Seurat
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(ggplot2)

####################
# Functions        #
####################

# read all lanes
read_all_lanes <- function(cellranger_lanes_dir, exclude_lanes = c(), min.cells = 3, min.features = 200){
  # start at null
  seurat_object <- NULL
  # get the subdirectories
  lanes <- list.dirs(cellranger_lanes_dir, recursive=F, full.names=F)
  # filter by exclusion lanes
  lanes <- setdiff(lanes, exclude_lanes)
  # grab all the data from each lane
  for(lane in lanes){
    seurat_object <- add_data(seurat_object, lane, cellranger_lanes_dir, min.cells, min.features)
  }
  # set the chemicality version of the lanes
  seurat_object$chem <- as.factor(ifelse(grepl(pattern = "^18", seurat_object$lane), "V2", "V3"))
  seurat_object$orig.ident <- as.factor(ifelse(grepl(pattern = "^18", seurat_object$lane), "stemi_v2", "stemi_v3"))
  return(seurat_object)
}

# add new data to Seurat object
add_data <- function(seurat_to_add_to = NULL, lane, base_counts_dir, min.cells = 3, min.features = 200) {
  print(lane)
  # get the counts
  counts_dir <- paste0(base_counts_dir, lane, "/outs/filtered_feature_bc_matrix/")
  # read the actual counts
  counts <- Read10X(counts_dir)
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # get the barcode without the '-X' appended, and the append the lane
  colnames(counts) <- paste0(substr(colnames(counts), 1, regexpr(last_dash_pos, colnames(counts))-1), "_", lane)
  # create some metadata, for now, we'll first just store the lane here
  metadata <- get_na_dataframe(c("batch","lane"),colnames(counts))
  metadata$lane = lane
  metadata$batch = lane
  # create the actual object
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = min.cells,
                                           min.features = min.features,
                                           project = "stemi",
                                           meta.data = metadata)
  # if we're starting from NULL, just return the current Seurat object
  if (is.null(seurat_to_add_to)){
    return(seurat_new)
  }
  # otherwise merge
  return(merge(seurat_to_add_to, seurat_new))
}


# plo the number of UMIs detected vs the mitochondrial percentage
plot_ncount_vs_mitopct <- function(seurat_object_metadata){
  # totally stolen from Harm
  p <- ggplot(seurat_object_metadata, aes(nCount_RNA, percent.mt)) + geom_hex(bins=100) +
    scale_fill_distiller(palette = "Spectral", name="Cell frequencies",
                         limits = c(0,100), oob = scales::squish) +
    ylab("Fraction mtDNA-encoded genes") + xlab("Number of UMIs") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "grey"),
          axis.text=element_text(size=12), axis.title=element_text(size=18))
  return(p)
}


# get a dataframe with NA values, with the given row and column names, needed if you want to add Seurat metadata, but don't have all values
get_na_dataframe <- function(colnames, rownames){
  # create matrix
  empty_matrix <- matrix(data = NA, nrow=length(rownames), ncol=length(colnames))
  # set the row and column names
  rownames(empty_matrix) <- rownames
  colnames(empty_matrix) <- colnames
  # convert to dataframe
  empty_frame <- as.data.frame(empty_matrix)
  return(empty_frame)
}


add_soup_assignments <- function(seurat_object, soup_dir, soup_append, batch_key='batch'){
  # grab the assignments
  soup_assignments <- get_soup_assignments(seurat_object, soup_dir, soup_append, batch_key)
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # add the lane/barcode combo
  rownames(soup_assignments) <- paste0(substr(soup_assignments$barcode, 1, regexpr(last_dash_pos, soup_assignments$barcode)-1), "_", soup_assignments$lane)
  # now add the information to the object
  for(column in colnames(soup_assignments)){
    #seurat_object <- AddMetaData(seurat_object, soup_assignments[column], paste("soup", column, sep = "_"))
    if(column != "lane" & column != "barcode")
      seurat_object <- AddMetaData(seurat_object, soup_assignments[column], column)
  }
  return(seurat_object)
}

# add scrublet doublet assignments
add_scrublet_assignments <- function(seurat_object, scrublet_loc){
  # read the output
  scrublet_output <- read.table(scrublet_loc, sep = '\t', header = T)
  # in case we accidently ran in twice on some samples
  scrublet_output <- unique(scrublet_output)
  rownames(scrublet_output) <- scrublet_output$lane_barcode
  # add the assignment
  seurat_object <- AddMetaData(seurat_object, scrublet_output['doublet'], 'scrublet_doublet')
  seurat_object <- AddMetaData(seurat_object, scrublet_output['doublet_score'], 'scrublet_dscore')
  return(seurat_object)
}

# add the assignments based on demuxlet demultiplexing
add_demux_assignments <- function(seurat_object, demux_dir, demux_append, batch_key='batch'){
  # grab the demux assignments
  demux_output_all <- get_demux_assignments(seurat_object, demux_dir, demux_append, batch_key)
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # add the lane/barcode combo
  rownames(demux_output_all) <- paste0(substr(demux_output_all$BARCODE, 1, regexpr(last_dash_pos, demux_output_all$BARCODE)-1), "_", demux_output_all$lane)
  # add the values
  seurat_object <- AddMetaData(seurat_object, demux_output_all['SNG.1ST'])
  seurat_object <- AddMetaData(seurat_object, demux_output_all['BEST'])
  seurat_object <- AddMetaData(seurat_object, demux_output_all['SNG.LLK1'])
  seurat_object <- AddMetaData(seurat_object, demux_output_all['LLK12'])
  return(seurat_object)
}


# grab a dataframe with the souporcell results for the lanes in the given object
get_soup_assignments <- function(seurat_object, soup_dir, soup_append, batch_key='batch'){
  # the method for getting demux output, also works for souporcell output
  soup_assignments <- get_demux_assignments(seurat_object, soup_dir, soup_append, batch_key)
  return(soup_assignments)
}

# grab a dataframe with the demux results for the lanes in the given object
get_demux_assignments <- function(seurat_object, demux_dir, demux_append, batch_key='batch'){
  # grab the unique lanes from the object
  unique_lanes <- unique(seurat_object@meta.data[[batch_key]])
  # we'll also create a frame to hold all the data
  demux_output_all <- NULL
  # go through those lanes
  for(i in 1:length(unique_lanes)){
    lane <- unique_lanes[i]
    print(lane)
    # create location of file
    demuxlet_output_file <- paste0(demux_dir, lane, demux_append)
    # read the actual file
    demuxlet_output <- read.table(demuxlet_output_file, header=T)
    # add the lane itself
    demuxlet_output$lane = lane
    # check if this was our first demux output
    if( is.null(demux_output_all)){
      demux_output_all <- demuxlet_output
    }
    # otherwise just append
    else{
      # add any missing columns in the new output (if less clusters than seen before)
      demuxlet_output[setdiff(colnames(demux_output_all), colnames(demuxlet_output))] <- NA
      # add any missing columns in the current output (if more clusters than seen before)
      demux_output_all[setdiff(colnames(demuxlet_output), colnames(demux_output_all))] <- NA
      demux_output_all <- rbind(demux_output_all, demuxlet_output)
    }
  }
  return(demux_output_all)
}

# do doublet detection and removal
remove_doublets <- function(seurat_object, detection_method="souporcell"){
  if(detection_method == "souporcell"){
    # the status keyword states whether this was a singlet or not
    seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data['status']!= "singlet", "possible_doublet")
    # then select by this new row
    seurat_object <- subset(seurat_object, subset = possible_doublet == F)
  }
  else if(detection_method == "demuxlet"){
    # include based on singlet likelyhood
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1 < 25)
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
  }
  return(seurat_object)
}

# add the stimulation tags, so UT, x hours after stim y, etc.
add_stim_tags <- function(seurat_object, stim_mapping_loc, assignment_key='assignment', batch_key='batch', tp_key='timepoint'){
  # add the column for the timepoints to the Seurat object
  seurat_object@meta.data[tp_key] <- NA
  # grab the timepoints
  tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
  # check the timepoints
  for(tp in colnames(tps)){
    print(tp)
    # check the lanes
    for(lane in rownames(tps)){
      # only apply if there are actually rows with lane in the object
      if(lane %in% seurat_object@meta.data[[batch_key]]){
        participants.as.string <- tps[lane,tp]
        # split the participant line by comma to get the participants for the timepoint
        participants.this.tp <- strsplit(participants.as.string, ",")
        # check if there are any cases with the combination of these participants with the timepoint (some participants had both v2 and v3 experiments)
        if(nrow(seurat_object@meta.data[seurat_object@meta.data[[assignment_key]] %in% unlist(participants.this.tp) 
                                        & seurat_object@meta.data[[batch_key]] == lane
                                        ,]) > 0){
          # set this timepoint for this lane combined with these participants
          seurat_object@meta.data[seurat_object@meta.data[[assignment_key]] %in% unlist(participants.this.tp) 
                                  & seurat_object@meta.data[[batch_key]] == lane
                                  ,][tp_key] <- tp
        }
      }
    }
  }
  return(seurat_object)
}


qc_table_to_pcts <- function(qc_numbers, lane_column='lane', reference='total') {
  with_pct <- apply(qc_numbers, 1, function(x) {
    # save the result as a list
    summaries <- list()
    # first the lane
    summaries[['lane']] <- x[[lane_column]]
    # get the total
    total <- as.numeric(x[[reference]])
    # get all the other columns
    others <- setdiff(names(x), c(lane_column))
    # check each column
    for (column in others) {
      # get the value
      value <- as.numeric(x[[column]])
      # calculate the percentage
      pct <- round(value / total, digits = 3) * 100
      # make into a string
      value_and_pct <- paste(value, ' ', '(', pct, '%)', sep = '')
      # put into result
      summaries[[column]] <- value_and_pct
    }
    return(data.frame(summaries))
  })
  return(do.call('rbind',with_pct))
}


####################
# main code        #
####################

# this is where our objects are on disk
object_loc <- "/groups/umcg-wijmenga/tmp02/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
# this is what we will save as
stemi_v3_raw_loc <- paste(object_loc, "stemi_v3_raw_samples.rds", sep = "/")
stemi_v3_filtered_loc <- paste(object_loc, "stemi_v3_filtered_samples.rds", sep = "/")
stemi_v2_raw_loc <- paste(object_loc, "stemi_v2_raw_samples.rds", sep = "/")
stemi_v2_filtered_loc <- paste(object_loc, "stemi_v2_filtered_samples.rds", sep = "/")
stemi_v3_raw_loc <- paste(object_loc, "stemi_v3_raw_samples_20201110.rds", sep = "/")
stemi_v3_filtered_loc <- paste(object_loc, "stemi_v3_filtered_samples_20201110.rds", sep = "/")
stemi_v2_raw_loc <- paste(object_loc, "stemi_v2_raw_samples_20201110.rds", sep = "/")
stemi_v2_filtered_loc <- paste(object_loc, "stemi_v2_filtered_samples_20201110.rds", sep = "/")
stemi_v2_normalized_loc <- paste(object_loc, "stemi_v2_normalized_samples_20201110.rds", sep = "/")
stemi_v3_normalized_loc <- paste(object_loc, "stemi_v3_normalized_samples_20201110.rds", sep = "/")

# this is where our cellranger outputs are
cellranger_lanes_dir <- "/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/alignment/hg19/cellranger_output/"

# read all the lanes
stemi <- read_all_lanes(cellranger_lanes_dir, exclude_lanes = c(), min.cells = 3, min.features = 200)
# separate the v2 and v3 ones
stemi_v2 <- subset(stemi, subset = chem == 'V2')
stemi_v3 <- subset(stemi, subset = chem == 'V3')
# save this raw file
saveRDS(stemi_v2, stemi_v2_raw_loc)
saveRDS(stemi_v3, stemi_v3_raw_loc)

# these are the soup pre- and appends
#soup_prepend <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/correlate_clusters/correlated_output/"
soup_prepend <- "/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/demultiplexing/souporcell_output/"
soup_append <- "_correlated.tsv"

# these are the demux pre- and appends
#demux_prepend <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/demuxlet/demuxlet_output/'
# demux_prepend <- '/groups/umcg-wijmenga/tmp02/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/demuxlet/demuxlet_output/'
# demux_append <- '_mmaf002.best'
# scrublet loc
#scrublet_v2_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/scrublet/scrublet_assignment_v2.tsv'
#scrublet_v3_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/scrublet/scrublet_assignment_v3.tsv'
# scrublet_v2_loc <- '/groups/umcg-wijmenga/tmp02/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/scrublet/scrublet_assignment_v2.tsv'
# scrublet_v3_loc <- '/groups/umcg-wijmenga/tmp02/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/scrublet/scrublet_assignment_v3.tsv'
# location of the simulation mapping
#stim_mapping_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/stemi-sampleIDs.txt'
stim_mapping_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/stemi-sampleIDs.txt'

# add souporcell assignments
stemi_v2 <- add_soup_assignments(stemi_v2, soup_prepend, soup_append)
# add demuxlet
# stemi_v2 <- add_demux_assignments(stemi_v2, demux_prepend, demux_append)
# add scrublet
# stemi_v2 <- add_scrublet_assignments(stemi_v2, scrublet_v2_loc)
# add stim tags
stemi_v2 <- add_stim_tags(stemi_v2, stim_mapping_loc = stim_mapping_loc, assignment_key = 'assignment_ll', tp_key = 'timepoint_ll')
# stemi_v2 <- add_stim_tags(stemi_v2, stim_mapping_loc = stim_mapping_loc, assignment_key = 'SNG.1ST', tp_key = 'timepoint.demux')

# save this raw file
saveRDS(stemi_v2, stemi_v2_raw_loc)

# calculate mt fraction
stemi_v2[["percent.mt"]] <- PercentageFeatureSet(stemi_v2, pattern = "^MT-")

# plot the mt fraction vs the gene count
mt_stemi_v2 <- plot_ncount_vs_mitopct(stemi_v2@meta.data) + ggtitle('total vs %mt gene count')
# saved as stemi_umi_vs_mtDNA_v2

# make some numbers for the QC
nrow(stemi_v2@meta.data)
# 66209
nrow(stemi_v2@meta.data[stemi_v2@meta.data$status == 'singlet', ]) # if it was a singlet
# 61435
nrow(stemi_v2@meta.data[
  stemi_v2@meta.data$status == 'singlet' & 
    stemi_v2@meta.data$nFeature_RNA > 200 &
    stemi_v2@meta.data$percent.mt < 8 &
    as.vector(stemi_v2@assays$RNA@counts['HBB', ]) < 10, ]
)
# 58418

# per lane as well
v2_unfiltered_per_lane <- data.frame(table(stemi_v2@meta.data[, 'lane']))
v2_singletfiltered_per_lane <- data.frame(table(stemi_v2@meta.data[stemi_v2@meta.data$status == 'singlet', 'lane']))
v2_singlet_and_qc_filtered_per_lane <- data.frame(table(stemi_v2@meta.data[stemi_v2@meta.data$status == 'singlet' & 
                                                                             stemi_v2@meta.data$nFeature_RNA > 200 &
                                                                             stemi_v2@meta.data$percent.mt < 8 &
                                                                             as.vector(stemi_v2@assays$RNA@counts['HBB', ]) < 10, 'lane']))
v2_lane_qc <- merge(
  merge(
    merge(v2_unfiltered_per_lane, v2_singletfiltered_per_lane, by = 'Var1'
    ), v2_singletfiltered_per_lane, by = 'Var1'
  ), v2_singlet_and_qc_filtered_per_lane, by = 'Var1')
colnames(v2_lane_qc) <- c('lane', 'cellranger', 'singlet filter', 'donor filter', 'rna filter')

# remove samples called as doublets by souporcell
stemi_v2 <- remove_doublets(stemi_v2, detection_method="souporcell")
# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
stemi_v2 <- subset(stemi_v2, subset = nFeature_RNA > 200 & percent.mt < 8 & HBB < 10)
# we're going to keep to the souporcell assignments for now
stemi_v2@meta.data$assignment.final <- stemi_v2@meta.data$assignment_ll
stemi_v2@meta.data$timepoint.final <- stemi_v2@meta.data$timepoint_ll

# save the preprocessed file
saveRDS(stemi_v2, stemi_v2_filtered_loc)

# add souporcell assignments
stemi_v3 <- add_soup_assignments(stemi_v3, soup_prepend, soup_append)
# fix for the missing GT
levels(stemi_v3@meta.data$assignment_ll) <- c(levels(stemi_v3@meta.data$assignment_ll), 'TEST_88')
stemi_v3@meta.data[!is.na(stemi_v3@meta.data$correlation_ll) & stemi_v3@meta.data$correlation_ll < 0.4, ]$assignment_ll <- 'TEST_88'
# add demuxlet
# stemi_v3 <- add_demux_assignments(stemi_v3, demux_prepend, demux_append)
# add scrublet
# stemi_v3 <- add_scrublet_assignments(stemi_v3, scrublet_v3_loc)
# add stim tags
stemi_v3 <- add_stim_tags(stemi_v3, stim_mapping_loc = stim_mapping_loc, assignment_key = 'assignment_ll', tp_key = 'timepoint_ll')
# stemi_v3 <- add_stim_tags(stemi_v3, stim_mapping_loc = stim_mapping_loc, assignment_key = 'SNG.1ST', tp_key = 'timepoint.demux')
# we're going to keep to the souporcell assignments for now
stemi_v3@meta.data$assignment.final <- stemi_v3@meta.data$assignment_ll
stemi_v3@meta.data$timepoint.final <- stemi_v3@meta.data$timepoint_ll

# save this raw file
saveRDS(stemi_v3, stemi_v3_raw_loc)

# calculate mt fraction
stemi_v3[["percent.mt"]] <- PercentageFeatureSet(stemi_v3, pattern = "^MT-")

# plot the mt fraction vs the gene count
mt_stemi_v3 <- plot_ncount_vs_mitopct(stemi_v3@meta.data) + ggtitle('total vs %mt gene count')
# saved as stemi_umi_vs_mtDNA_v3 (10x10)

# make some numbers for the QC
nrow(stemi_v3@meta.data)
# 71881
nrow(stemi_v3@meta.data[stemi_v3@meta.data$status == 'singlet', ]) # if it was a singlet
# 62931
nrow(stemi_v3@meta.data[
  stemi_v3@meta.data$status == 'singlet' &
    stemi_v3@meta.data$assignment.final != 'TEST_60' & 
    stemi_v3@meta.data$assignment.final != 'TEST_62' & 
    !(stemi_v3@meta.data$assignment.final == 'TEST_68' & stemi_v3@meta.data$timepoint.final == 't24h'), ]
)
# 55343
nrow(stemi_v3@meta.data[
  stemi_v3@meta.data$status == 'singlet' &
    stemi_v3@meta.data$assignment.final != 'TEST_60' & 
    stemi_v3@meta.data$assignment.final != 'TEST_62' & 
    !(stemi_v3@meta.data$assignment.final == 'TEST_68' & stemi_v3@meta.data$timepoint.final == 't24h') & 
    stemi_v3@meta.data$nFeature_RNA > 200 &
    stemi_v3@meta.data$percent.mt < 15 &
    as.vector(stemi_v3@assays$RNA@counts['HBB', ]) < 10, ]
)
# 37577

# per lane as well
v3_unfiltered_per_lane <- data.frame(table(stemi_v3@meta.data[, 'lane']))
v3_singletfiltered_per_lane <- data.frame(table(stemi_v3@meta.data[stemi_v3@meta.data$status == 'singlet', 'lane']))
v3_singlet_and_included_per_lane <- data.frame(table(stemi_v3@meta.data[
  stemi_v3@meta.data$status == 'singlet' &
    stemi_v3@meta.data$assignment.final != 'TEST_60' & 
    stemi_v3@meta.data$assignment.final != 'TEST_62' & 
    !(stemi_v3@meta.data$assignment.final == 'TEST_68' & stemi_v3@meta.data$timepoint.final == 't24h'), 'lane']))
v3_singlet_included_and_qc_filtered_per_lane <- data.frame(table(stemi_v3@meta.data[
  stemi_v3@meta.data$status == 'singlet' &
    stemi_v3@meta.data$assignment.final != 'TEST_60' & 
    stemi_v3@meta.data$assignment.final != 'TEST_62' & 
    !(stemi_v3@meta.data$assignment.final == 'TEST_68' & stemi_v3@meta.data$timepoint.final == 't24h') & 
    stemi_v3@meta.data$nFeature_RNA > 200 &
    stemi_v3@meta.data$percent.mt < 15 &
    as.vector(stemi_v3@assays$RNA@counts['HBB', ]) < 10, 'lane']))
v3_lane_qc <- merge(
  merge(
    merge(
      v3_unfiltered_per_lane, v3_singletfiltered_per_lane, by = 'Var1'
    ), v3_singlet_and_included_per_lane, by = 'Var1'
  ), v3_singlet_included_and_qc_filtered_per_lane, by = 'Var1'
)
colnames(v3_lane_qc) <- c('lane', 'cellranger', 'singlet filter', 'donor filter', 'rna filter')

# remove samples called as doublets by souporcell
stemi_v3 <- remove_doublets(stemi_v3, detection_method="souporcell")

# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
stemi_v3 <- subset(stemi_v3, subset = nFeature_RNA > 200 & percent.mt < 15 & HBB < 10)
# save the preprocessed file
saveRDS(stemi_v3, stemi_v3_filtered_loc)

# remove some samples that should not be in there
stemi_v3 <- subset(stemi_v3, subset = assignment.final != 'TEST_60' & assignment.final != 'TEST_62')
stemi_v3 <- stemi_v3[, !(stemi_v3@meta.data$assignment.final == 'TEST_68' & stemi_v3@meta.data$timepoint.final == 't24h')]

# save the QC
stemi_qc_table_lanes <- qc_table_to_pcts(rbind(v2_lane_qc, v3_lane_qc), reference = 'cellranger')

# do normalization
stemi_v2 <- NormalizeData(stemi_v2)
stemi_v2 <- SCTransform(stemi_v2, vars.to.regress = c('percent.mt'))
stemi_v3 <- NormalizeData(stemi_v3)
stemi_v3 <- SCTransform(stemi_v3, vars.to.regress = c('percent.mt'))
saveRDS(stemi_v2, stemi_v2_normalized_loc)
saveRDS(stemi_v3, stemi_v3_normalized_loc)