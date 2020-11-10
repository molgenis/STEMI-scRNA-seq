####################
# libraries        #
####################

library(Seurat)

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
  return(seurat_object)
}

# add new data to Seurat object
add_data <- function(seurat_to_add_to = NULL, lane, base_counts_dir, min.cells = 3, min.features = 200) {
  print(lane)
  # get the counts
  counts_dir <- paste0(base_counts_dir, lane, "/outs/filtered_feature_bc_matrix/")
  # read the actual counts
  counts <- Read10X(counts_dir)
  
  # to keep the barcodes unique, we're appending a number, lets keep increasing that number it not the first one
  barcode_append = 1
  if (is.null(seurat_to_add_to) == F){
    barcode_append = length(unique(seurat_to_add_to@meta.data$batch)) + 1
  }
  # append that number
  colnames(counts) <- paste0(colnames(counts), "-1-",barcode_append)
  
  # create some metadata, for now, we'll first just store the lane here
  metadata <- get_na_dataframe(c("batch","lane"),colnames(counts))
  metadata$lane = lane
  metadata$batch = lane
  # create the actual object
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = min.cells,
                                           min.features = min.features,
                                           project = "1M_cells",
                                           meta.data = metadata)
  # if we're starting from NULL, just return the current Seurat object
  if (is.null(seurat_to_add_to)){
    return(seurat_new)
  }
  # otherwise merge
  return(merge(seurat_to_add_to, seurat_new))
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

# add the lanes back to the object by reading a tsv with the barcodes as they were in the h5ad, and the lane beloninging to it
readd_lanes <- function(seurat_object, lane_loc){
  # read the lanes beloninging to the barcodes back (file was made from h5ad file)
  lanes <- read.table(lane_loc, header = T, row.names = 1)
  # add the information back to the object
  AddMetaData(seurat_object, lanes, col.name="batch")
}

# add the experiment tags, the experiment number beloninging to an LLID
add_exp_tags <- function(seurat_object, exp_to_ll_loc, assignment_key='assignment'){
  # add the column for expIDs to the Seurat object
  seurat_object@meta.data$exp.id <- NA
  # grab the mapping of expnr to ll
  exp_mapping = read.table(exp_to_ll_loc, header = T, stringsAsFactors = F)
  # check for each LL ID what the experiment number was
  for(lld in exp_mapping$LLD.ID){
    # grab the experiment number
    exp_nr <- exp_mapping[exp_mapping$LLD.ID == lld,]$ExpNr
    # no need to overwrite if the value was already NA
    if(is.na(exp_nr) == F & nrow(seurat_object@meta.data[seurat_object@meta.data[assignment_key] == lld & is.na(seurat_object@meta.data[assignment_key])==F,]) > 0){
      # set the experiment number for all those LL IDs
      seurat_object@meta.data[seurat_object@meta.data[assignment_key] == lld & is.na(seurat_object@meta.data[assignment_key])==F,]$exp.id <- exp_nr
    }
  }
  return(seurat_object)
}

# add the stimulation tags, so UT, x hours after stim y, etc.
add_stim_tags <- function(seurat_object, stim_mapping_loc, assignment_key='exp.id', batch_key='batch', tp_key='timepoint'){
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

# add base barcodes to seurat metadata, to without the last append that scanpy added
add_base_barcodes <- function(seurat_object, colname="barcode_base"){
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # get the full barcodes
  barcodes_full <- rownames(seurat_object@meta.data)
  # create the cut barcodes
  barcodes_cut <- substr(barcodes_full, 1, regexpr(last_dash_pos, barcodes_full)-1)
  # add this to the Seurat Object
  seurat_object@meta.data[colname] <- barcodes_cut
  return(seurat_object)
}

# add the base barcode + lane as value in seurat metadata (requires for some matching stuff)
add_lane_barcode_combo <- function(seurat_object, colname="barcode_lane"){
  # add the base barcode if it was not already there
  if("barcode_base" %in% colnames(seurat_object@meta.data) == F){
    seurat_object <- add_base_barcodes(seurat_object)
  }
  # combine the lane and barcode into a string
  seurat_object@meta.data[colname] <- paste0(seurat_object@meta.data$batch,"_",seurat_object$barcode_base)
  return(seurat_object)
}

# add the bare barcode, so without a number append
add_bare_barcodes <- function(seurat_object, colname="barcode_bare"){
  # create a regex to get the first index of -
  last_dash_pos <- "\\-"
  # get the full barcodes
  barcodes_full <- rownames(seurat_object@meta.data)
  # create the cut barcodes
  barcodes_cut <- substr(barcodes_full, 1, regexpr(last_dash_pos, barcodes_full)-1)
  # add this to the Seurat Object
  seurat_object@meta.data[colname] <- barcodes_cut
  return(seurat_object)
}

# add a combination of the lane with the bare barcode e.g. ACTGGACA_170326
add_lane_bare_barcode_combo <- function(seurat_object, colname="bare_barcode_lane"){
  if("barcode_bare" %in% colnames(seurat_object@meta.data) == F){
    seurat_object <- add_bare_barcodes(seurat_object)
  }
  # combine the lane and barcode into a string
  seurat_object@meta.data[colname] <- paste0(seurat_object$barcode_bare,"_",seurat_object@meta.data$batch)
  return(seurat_object)
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

# where to store the objects
object_loc <- "/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
# this directory contains the cellranger outputs
cellranger_lanes_dir <- "/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/alignment/cellranger_output/"
# we'll save some plots here
plot_dir = "/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/"
# excluding some lanes
exclude_lanes <- c("181010_lane3", "181011_lane3", "181017_lane3", "181105_lane3", "181106_lane3", "181121_lane3", "181122_lane3", "181213_lane4", "181214_lane3","181024_lane1","181024_lane2","181024_lane3","190101_lane1","190101_lane2")
# include HCs
include_hc_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/metadata/hc-sampleIDs.txt'

# this directory houses the barcode assignments to the participants based on souporcell
base_soup_dir <- "/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/correlate_clusters/correlated_output/"
# this is the append of the barcode assignments souporcell file
soup_extension <- "_correlated.tsv"
# this directory houses the barcode assignments to the participants based on demuxlet
base_demux_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/demux_output_cytosnp/"
#this is the append of the barcode assignments demux file
demux_extension <- "_sorted_hfixed.best"

# this file contains the lanes as rownames, timepoints as colnames and the partipants in the cells, split by a comma
stim_mapping_loc <- "/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_lane_to_tp.tsv"
# this file contains the expnr to ll id
exp_to_ll_loc <- "/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_to_ll.tsv"


# read all the lanes
cells_1M <- read_all_lanes(cellranger_lanes_dir, exclude_lanes = exclude_lanes, min.cells = 3, min.features = 200)
# set better colnames
cells_1M <- add_bare_barcodes(cells_1M)
cells_1M <- add_lane_bare_barcode_combo(cells_1M)
cells_1M <- RenameCells(cells_1M, new.names = cells_1M@meta.data$bare_barcode_lane)
# add the identities
#cells_1M <- add_demux_assignments(cells_1M, base_demux_dir, demux_extension)
# add the exp nr
#cells_1M <- add_exp_tags(cells_1M, exp_to_ll_loc, assignment_key='SNG.1ST')
# we want to add the exp tag based on souporcell as well, so let's store the exp column under an other name as well
#cells_1M@meta.data$exp.id.demux <- cells_1M@meta.data$exp.id
# add the condition
#cells_1M <- add_stim_tags(cells_1M, stim_mapping_loc)
# we want to add the condition tag based on souporcell as well, so let's store the exp column under an other name as well
#cells_1M@meta.data$timepoint.demux <- cells_1M@meta.data$timepoint
# add souporcell assignments
cells_1M <- add_soup_assignments(cells_1M, base_soup_dir, soup_extension)
# add the exp nr for souporcell this time
cells_1M <- add_exp_tags(cells_1M, exp_to_ll_loc, assignment_key='assignment_ll')
cells_1M@meta.data$exp.id.ll <- cells_1M@meta.data$exp.id
# add the condition again, but now the exp nr is the ll one
cells_1M <- add_stim_tags(cells_1M, stim_mapping_loc)
cells_1M@meta.data$timepoint.ll <- cells_1M@meta.data$timepoint
# remove the doublets
cells_1M <- remove_doublets(cells_1M)
# set the final assignments
cells_1M@meta.data$assignment.final <- cells_1M@meta.data$assignment.ll
cells_1M@meta.data$timepoint.final <- cells_1M@meta.data$timepoint.ll
cells_1M@meta.data$exp.id.final <- cells_1M@meta.data$exp.id.ll
# grab just UT
cells_1M <- subset(cells_1M, subset = timepoint.final == 'UT')
# grab just the participants that are HC
include_hc_list <- read.table(include_hc_loc)
cells_1M <- subset(cells_1M, subset = assignment.final %in% include_hc_list$V1)
# split on chem
HC_v2 <- subset(cells_1M, subset = chem == 'V2')
HC_v3 <- subset(cells_1M, subset = chem == 'V3')
# calculate mt fraction
HC_v2[["percent.mt"]] <- PercentageFeatureSet(HC_v2, pattern = "^MT-")
HC_v3[["percent.mt"]] <- PercentageFeatureSet(HC_v3, pattern = "^MT-")
# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
HC_v2 <- subset(HC_v2, subset = nFeature_RNA > 200 & percent.mt < 8 & HBB < 10)
HC_v3 <- subset(HC_v3, subset = nFeature_RNA > 200 & percent.mt < 15 & HBB < 10)
# do normalization
HC_v2 <- SCTransform(HC_v2, vars.to.regress = c('percent.mt'))
HC_v3 <- SCTransform(HC_v3, vars.to.regress = c('percent.mt'))
# save the HC objects
saveRDS(HC_v2, paste(object_loc, 'HC_v2_20201110.rds', sep = ''))
saveRDS(HC_v3, paste(object_loc, 'HC_v3_20201110.rds', sep = ''))
