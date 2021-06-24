#
# header describing file
#


# load the libraries
library(Seurat)
library(Matrix)
library(CellChat)
library(patchwork)

# set some options
options(stringsAsFactors = FALSE)

# functions
init_cellchat_object <- function(seurat_object, assay = 'SCT', slot = 'data', ident = 'cell_type_lowerres'){
  # set the default assay
  DefaultAssay(seurat_object) <- assay
  # extract the data
  data.input <- GetAssayData(seurat_object, assay = assay, slot = slot)
  meta <- seurat_object@meta.data
  # create the object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = ident)
  return(cellchat)
}

preprocess_cellchat_object <- function(chat_object, nthreads=8){
  # set the database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  chat_object@DB <- CellChatDB.use
  chat_object <- subsetData(chat_object) # This step is necessary even if using the whole database
  # set multithreading options
  future::plan("multiprocess", workers = nthreads) # do parallel
  # get genes
  chat_object <- identifyOverExpressedGenes(chat_object)
  chat_object <- identifyOverExpressedInteractions(chat_object)
  # project gene expression data onto PPI network (optional)
  chat_object <- projectData(chat_object, PPI.human)
  return(chat_object)
}

inference_communication_network <- function(chat_object, min.cells = 10){
  # Compute the communication probability and infer cellular communication network
  chat_object <- computeCommunProb(chat_object)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  chat_object <- filterCommunication(chat_object, min.cells = min.cells)
  # Infer the cell-cell communication at a signaling pathway level
  chat_object <- computeCommunProbPathway(chat_object)
  # Calculate the aggregated cell-cell communication network
  chat_object <- aggregateNet(chat_object)
  return(chat_object)
}

plot_communication_network <- function(chat_object){
  groupSize <- as.numeric(table(chat_object@idents))
  par(mfrow = c(1,2), xpd=TRUE)
  netVisual_circle(chat_object@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
  netVisual_circle(chat_object@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
}

plot_communication_network_separate <- function(chat_object, slot = 'weight'){
  # grab the number of groups
  groupSize <- as.numeric(table(chat_object@idents))
  # grab the right slot
  mat <- chat_object@net[[slot]]
  # make a perfect square of the plot
  par(mfrow = c(ceiling(sqrt(nrow(mat))),ceiling(sqrt(nrow(mat)))), xpd=TRUE)
  # make the plots
  for (i in 1:nrow(mat)) {
    mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    mat2[i, ] <- mat[i, ]
    netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
}


# set paths
object_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
cardio.integrated.loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')

# read object
cardio.integrated <- readRDS(cardio.integrated.loc)
# subset to specific timepoint and timepoint
cardio.stemi_v2.baseline <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v2' & cardio.integrated@meta.data$timepoint.final == 'Baseline']
# create cell type
cardio.stemi_v2.baseline.chat <- init_cellchat_object(cardio.stemi_v2.baseline)
cardio.stemi_v2.baseline.chat <- preprocess_cellchat_object(cardio.stemi_v2.baseline.chat)
cardio.stemi_v2.baseline.chat <- inference_communication_network(cardio.stemi_v2.baseline.chat)
# plot the interactions
plot_communication_network_separate(cardio.stemi_v2.baseline.chat)

