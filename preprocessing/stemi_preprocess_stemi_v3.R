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
                                           project = "cardio_v3",
                                           meta.data = metadata)
  # if we're starting from NULL, just return the current Seurat object
  if (is.null(seurat_to_add_to)){
    return(seurat_new)
  }
  # otherwise merge
  return(merge(seurat_to_add_to, seurat_new))
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



####################
# main code        #
####################

# this is where our objects are on disk
object_loc <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
# this is what we will save as
stemi_v3_raw_loc <- paste(object_loc, "stemi_v3_raw_samples.rds", sep = "/")
stemi_v3_filtered_loc <- paste(object_loc, "stemi_v3_filtered_samples.rds", sep = "/")

# this is where our cellranger outputs are
cellranger_lanes_dir_v3 <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cellranger_output/"

# read all the lanes
stemi_v3 <- read_all_lanes(cellranger_lanes_dir_v3, exclude_lanes = exclude_lanes, min.cells = 3, min.features = 200)
# save this raw file
saveRDS(stemi_v3, stemi_v3_raw_loc)

# these are the soup pre- and appends
soup_prepend <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/soupor_doublet/doublet_output/"
soup_append <- "_unremapped_doublets.tsv"

# add (preliminary) souporcell assignments
stemi_v3 <- add_soup_assignments(stemi_v3, soup_prepend, soup_append)

# remove samples called as doublets by souporcell
stemi_v3 <- remove_doublets(stemi_v3, detection_method="souporcell")

# calculate mt fraction
stemi_v3[["percent.mt"]] <- PercentageFeatureSet(stemi_v3, pattern = "^MT-")
# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
stemi_v3 <- subset(stemi_v3, subset = nFeature_RNA > 200 & percent.mt < 15 & HBB < 10)
# save the preprocessed file
saveRDS(stemi_v3, stemi_v3_filtered_loc)
