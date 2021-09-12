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

inference_communication_network <- function(chat_object, min.cells=10, thresh=1){
  # Compute the communication probability and infer cellular communication network
  chat_object <- computeCommunProb(chat_object)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  chat_object <- filterCommunication(chat_object, min.cells = min.cells)
  # Infer the cell-cell communication at a signaling pathway level
  chat_object <- computeCommunProbPathway(chat_object, thresh = thresh)
  # Calculate the aggregated cell-cell communication network
  chat_object <- aggregateNet(chat_object, thresh = thresh)
  return(chat_object)
}

do_default_cellchat_workflow <- function(seurat_object, assay = 'SCT', slot = 'data', ident = 'cell_type_lowerres', min.cells=10, nthreads=8, thresh=0.05){
  # go through the steps
  chat_object <- init_cellchat_object(seurat_object, assay = assay, slot = slot, ident = ident)
  chat_object <- preprocess_cellchat_object(chat_object, nthreads = nthreads)
  chat_object <- inference_communication_network(chat_object, min.cells = min.cells, thresh = thresh)
  return(chat_object)
}

do_default_cellchat_workflow_per_timepoint <- function(seurat_object, timepoint.column='timepoint.final', assay = 'SCT', slot = 'data', ident = 'cell_type_lowerres', min.cells=10, nthreads=8, thresh=0.05){
  # save the objects in a list
  chat_per_timepoint <- list()
  # check each cell type
  for(timepoint in unique(seurat_object@meta.data[[timepoint.column]])){
    # subset to that timepoint
    seurat_timepoint <- seurat_object[, seurat_object@meta.data[[timepoint.column]] == timepoint]
    # go through the work flow
    chat_timepoint <- do_default_cellchat_workflow(seurat_timepoint, assay = assay, slot = slot, ident = ident, min.cells=min.cells, nthreads=nthreads, thresh = thresh)
    # add to list
    chat_per_timepoint[[timepoint]] <- chat_timepoint
  }
  return(chat_per_timepoint)
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

plot_all_communications_networks_separate <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10, slot = 'weight'){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # plot each chem
    for(chem in chems){
      # paste the full output loc
      output_loc_full <- paste(output_loc, 'communications_separate_', timepoint, '_', chem, '_', slot, '.', file_type, sep = '')
      print(output_loc_full)
      # init where we will save
      if(file_type == 'pdf'){
        pdf(output_loc_full, width = width, height = height)
      }
      else if(file_type == 'png'){
        png(output_loc_full, width = width, height = height)
      }
      else{
        print('unknown file type, doing pdf instead')
        pdf(output_loc_full, width = width, height = height)
      }
      try({
      # grab the number of groups
      groupSize <- as.numeric(table(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@idents))
      # grab the right slot
      mat <- chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net[[slot]]
      # make a perfect square of the plot
      par(mfrow = c(ceiling(sqrt(nrow(mat))),ceiling(sqrt(nrow(mat)))), xpd=TRUE)
      # make the plots
      for (i in 1:nrow(mat)) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = paste(rownames(mat)[i], 'in', chem, timepoint))
      }
      })
      dev.off()
    }
  }
  par(mfrow=c(1,1))
}

plot_all_communications_networks <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # paste the full output loc
    output_loc_full <- paste(output_loc, 'communications_', timepoint, '.', file_type, sep = '')
    print(output_loc_full)
    # init where we will save
    if(file_type == 'pdf'){
      pdf(output_loc_full, width = width, height = height)
    }
    else if(file_type == 'png'){
      png(output_loc_full, width = width, height = height)
    }
    else{
      print('unknown file type, doing pdf instead')
      pdf(output_loc_full, width = width, height = height)
    }
    try({
      # plot each chem
      for(chem in chems){
        groupSize <- as.numeric(table(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@idents))
        netVisual_circle(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Number of interactions", 'in', chem, timepoint))
        netVisual_circle(chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = paste("Interaction weights/strength", 'in', chem, timepoint))
      }
    })
    dev.off()
  }
}

plot_all_aggregate_pathways <- function(chat_object_per_timepoint_and_chem, output_loc='', file_type='pdf', chems=c('V2', 'V3'), width=10, height=10){
  # plot each timepoint
  for(timepoint in names(chat_object_per_timepoint_and_chem[[chems[[1]]]])){
    # plot each chem
    for(chem in chems){
      # paste the full output loc
      output_loc_full <- paste(output_loc, 'pathways_', timepoint, '_', chem, '.', file_type, sep = '')
      print(output_loc_full)
      # init where we will save
      if(file_type == 'pdf'){
        pdf(output_loc_full, width = width, height = height)
      }
      else if(file_type == 'png'){
        png(output_loc_full, width = width, height = height)
      }
      else{
        print('unknown file type, doing pdf instead')
        pdf(output_loc_full, width = width, height = height)
      }
      try({
        par(mfrow=c(3,2))
        for(pathway in chat_object_per_timepoint_and_chem[[chem]][[timepoint]]@netP$pathways){
          netVisual_aggregate(chat_object_per_timepoint_and_chem[[chem]][[timepoint]], signaling = c(pathway), layout = "circle")
        }
      })
      dev.off()
    }
  }
}

compare_conditions <- function(all_conditions_list, condition.1, condition.2){
  # subset to only the conditions we want
  conditions.list <- list(condition.1 = all_conditions_list[[condition.1]], condition.2 = all_conditions_list[[condition.2]])
  conditions <- mergeCellChat(conditions.list, add.names = names(conditions.list))
  # compare interactions and strength numbers
  gg1 <- compareInteractions(conditions, show.legend = F, group = c(1,2))
  gg2 <- compareInteractions(conditions, show.legend = F, group = c(1,2), measure = "weight")
  # check strengths and weights as circle plots
  gg3 <- netVisual_diffInteraction(conditions, weight.scale = T)
  gg4 <- netVisual_diffInteraction(conditions, weight.scale = T, measure = "weight")
  # check strengths and weights as heatmaps
  gg5 <- netVisual_heatmap(conditions)
  gg6 <- netVisual_heatmap(conditions, measure = "weight")
  # compute net neutrality score to show differences
  conditions.list[[1]] <- netAnalysis_computeCentrality(conditions.list[[1]])
  conditions.list[[2]] <- netAnalysis_computeCentrality(conditions.list[[2]])
  # visualize in 2d space, #first get the number of links in both sets
  conditions.list.num.link <- sapply(conditions.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
  # set a max and min, so the plots will have the same scale
  conditions.list.weight.MinMax <- c(min(conditions.list.num.link), max(conditions.list.num.link))
  # make the plots
  gg7 <- netAnalysis_signalingRole_scatter(conditions.list[[1]], title = names(conditions.list)[1], weight.MinMax = conditions.list.weight.MinMax)
  gg8 <- netAnalysis_signalingRole_scatter(conditions.list[[2]], title = names(conditions.list)[2], weight.MinMax = conditions.list.weight.MinMax)
  # put everything in a nice list
  all_plots <- list('interactions_number' = gg1,
                    'interactions_weight' = gg2,
                    'circle_number' = gg3,
                    'circle_weight' = gg4,
                    'heatmap_number' = gg5,
                    'heatmap_weight' = gg6,
                    'scatter_cond1' = gg7,
                    'scatter_cond2' = gg8)
  return(all_plots)
}



# set paths
object_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/'
cardio.integrated.loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')
chat.all.loc <- paste(object_loc, 'cardio.chat.all.20210301.rds', sep = '')

# read object
cardio.integrated <- readRDS(cardio.integrated.loc)
# use only our major cell types
cardio.integrated <- cardio.integrated[ , cardio.integrated@meta.data$cell_type_lowerres %in% c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')]
cardio.integrated@meta.data$cell_type_lowerres <- droplevels(cardio.integrated@meta.data$cell_type_lowerres)

# subset to specific timepoint and timepoint
cardio.stemi_v2.baseline <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v2' & cardio.integrated@meta.data$timepoint.final == 'Baseline']
# create cell type
cardio.stemi_v2.baseline.chat <- init_cellchat_object(cardio.stemi_v2.baseline)
cardio.stemi_v2.baseline.chat <- preprocess_cellchat_object(cardio.stemi_v2.baseline.chat)
cardio.stemi_v2.baseline.chat <- inference_communication_network(cardio.stemi_v2.baseline.chat)
# plot the interactions
plot_communication_network_separate(cardio.stemi_v2.baseline.chat)

# same for v3
# subset to specific timepoint and timepoint
cardio.stemi_v3.baseline <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v3' & cardio.integrated@meta.data$timepoint.final == 'Baseline']
# create cell type
cardio.stemi_v3.baseline.chat <- init_cellchat_object(cardio.stemi_v3.baseline)
cardio.stemi_v3.baseline.chat <- preprocess_cellchat_object(cardio.stemi_v3.baseline.chat)
cardio.stemi_v3.baseline.chat <- inference_communication_network(cardio.stemi_v3.baseline.chat)
# now just for each timepoint
cardio.chem_v2.all_list <- do_default_cellchat_workflow_per_timepoint(cardio.integrated[, cardio.integrated@meta.data$chem == 'V2'], nthreads = 1, thresh = 1)
cardio.chem_v3.all_list <- do_default_cellchat_workflow_per_timepoint(cardio.integrated[, cardio.integrated@meta.data$chem == 'V3'], nthreads = 1, thresh = 1)
# save in list
cardio.list_per_chem <- list('V2' = cardio.chem_v2.all_list, 'V3' = cardio.chem_v3.all_list)
saveRDS(cardio.list_per_chem, chat.all.loc)
plot_all_communications_networks_separate(cardio.list_per_chem, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/cellchat/')
plot_all_aggregate_pathways(cardio.list_per_chem, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_cell_interactions/plots/cellchat/')


# start comparing the conditions
cardio.stemi_v3.t24h <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v3' & cardio.integrated@meta.data$timepoint.final == 't24h']
# create cell type
cardio.stemi_v3.t24h.chat <- init_cellchat_object(cardio.stemi_v3.t24h)
cardio.stemi_v3.t24h.chat <- preprocess_cellchat_object(cardio.stemi_v3.t24h.chat)
cardio.stemi_v3.t24h.chat <- inference_communication_network(cardio.stemi_v3.t24h.chat)
# and 8w
cardio.stemi_v3.t8w <- cardio.integrated[, cardio.integrated@meta.data$orig.ident == 'stemi_v3' & cardio.integrated@meta.data$timepoint.final == 't8w']
cardio.stemi_v3.t8w.chat <- init_cellchat_object(cardio.stemi_v3.t8w)
cardio.stemi_v3.t8w.chat <- preprocess_cellchat_object(cardio.stemi_v3.t8w.chat)
cardio.stemi_v3.t8w.chat <- inference_communication_network(cardio.stemi_v3.t8w.chat)
# finally let's do UT
cardio.stemi_v3.UT <- cardio.integrated[, cardio.integrated@meta.data$chem == 'V3' & cardio.integrated@meta.data$timepoint.final == 'UT']
cardio.stemi_v3.UT.chat <- init_cellchat_object(cardio.stemi_v3.UT)
cardio.stemi_v3.UT.chat <- preprocess_cellchat_object(cardio.stemi_v3.UT.chat)
cardio.stemi_v3.UT.chat <- inference_communication_network(cardio.stemi_v3.UT.chat)
# put everything in the list
cardio.stemi_v3.all.list <- list('Baseline' = cardio.stemi_v3.baseline.chat, 't24h' = cardio.stemi_v3.t24h.chat, 't8w' = cardio.stemi_v3.t8w.chat, 'UT' = cardio.stemi_v3.UT.chat)
# do analyses we care about
cardio.stemi.v3.chatresult.Baseline_t24h <- compare_conditions(cardio.stemi_v3.all.list, 'Baseline', 't24h')
cardio.stemi.v3.chatresult.Baseline_t8w <- compare_conditions(cardio.stemi_v3.all.list, 'Baseline', 't8w')
cardio.stemi.v3.chatresult.UT_Baseline <- compare_conditions(cardio.stemi_v3.all.list, 'UT', 'Baseline')


