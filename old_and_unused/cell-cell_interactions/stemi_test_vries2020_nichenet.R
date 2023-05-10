library(Seurat)
library(SeuratDisk)
library(nichenetr)
library(tidyverse)


# update to newer version to get more usable 
vries2020 <- readRDS('/data/scRNA/Seurat/objects/pilot4.Rds')
vries2020 <- UpdateSeuratObject(vries2020)
# try to turn ensemble IDs into gene symbols
symbols.to.ensg.mapping <- "/data/cardiology/eQTL_mapping/features_v3.tsv"
genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
genes$V2 <- gsub("_", "-", make.unique(genes$V2))
# grab the counts matrix
vries2020_counts <- as.matrix(vries2020@assays$RNA@counts)
symbols_rows <- genes[match(rownames(vries2020_counts), genes$V1),"V2"]
rownames(vries2020_counts) <- symbols_rows
vries2020 <- CreateSeuratObject(vries2020_counts, project='pp2020', meta.data=vries2020@meta.data)

# normalize pilot3 in the same ways as the reference was, with SCTransform
vries2020 <- SCTransform(vries2020)
vries2020 <- RunPCA(vries2020)

# read the reference cell typing object
reference <- LoadH5Seurat('/data/scRNA/Seurat/objects/pbmc_multimodal.h5seurat')

# find anchors between the reference and vries2020
anchors <- FindTransferAnchors(
  reference = reference,
  query = vries2020,
  normalization.method = "SCT",
  reference.reduction = "spca",
  dims = 1:50,
  recompute.residuals = FALSE # workaround for Seurat bug
)

# project the query data onto the reference structure
vries2020 <- MapQuery(
  anchorset = anchors,
  query = vries2020,
  reference = reference,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    predicted_ADT = "ADT"
  ),
  reference.reduction = "spca", 
  reduction.model = "wnn.umap"
)

# see what manual clustering tells us
vries2020 <- FindNeighbors(vries2020, dims = 1:30, verbose = TRUE)
vries2020 <- FindClusters(vries2020, resolution = 1.2, verbose = TRUE)
vries2020 <- RunUMAP(vries2020, dims = 1:30, verbose = TRUE)

# we will use the l2 predicted cell types (for the most part), for the classification of the object
vries2020@meta.data$cell_type <- as.character(vries2020@meta.data$predicted.celltype.l2)
# add another string to denote these as unclassified
vries2020@meta.data[is.na(vries2020@meta.data$cell_type), 'cell_type'] <- 'unclassified'

# harmonize the cell types again
vries2020@meta.data[vries2020@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
vries2020@meta.data[vries2020@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
vries2020@meta.data[vries2020@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
vries2020@meta.data[vries2020@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
vries2020@meta.data[vries2020@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
vries2020@meta.data[vries2020@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
vries2020@meta.data[vries2020@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
vries2020@meta.data[vries2020@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
# spaces in variables is inconvenient in a lot of places, so we'll replace these with underscores
vries2020@meta.data$cell_type <- gsub(' ', '_', vries2020@meta.data$cell_type)

# we will define a lower resolution cell type as well, we need to create some groupings for this
cd4t <- c('Treg', 'CD4_Naive', 'CD4_TCM', 'CD4_TEM', 'CD4_CTL', 'CD4_Proliferating')
cd8t <- c('MAIT', 'CD8_Naive', 'CD8_TCM', 'CD8_TEM', 'CD8_Proliferating')
t_other <- c('dnT', 'gdT', 'ILC')
nk <- c('NKdim', 'NKbright', 'NK_Proliferating')
monocyte <- c('cMono', 'ncMono')
dc <- c('cDC1', 'cDC2', 'pDC', 'ASDC')
b <- c('B_naive', 'B_intermediate', 'B_memory')
# add the new column by copying the higher res first, 
vries2020@meta.data$cell_type_lowerres <- vries2020@meta.data$cell_type
# in this new column, overwrite them to have the lower resolution
vries2020@meta.data[vries2020@meta.data$cell_type %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
vries2020@meta.data[vries2020@meta.data$cell_type %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
vries2020@meta.data[vries2020@meta.data$cell_type %in% t_other, ]$cell_type_lowerres <- 'T_other'
vries2020@meta.data[vries2020@meta.data$cell_type %in% nk, ]$cell_type_lowerres <- 'NK'
vries2020@meta.data[vries2020@meta.data$cell_type %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
vries2020@meta.data[vries2020@meta.data$cell_type %in% dc, ]$cell_type_lowerres <- 'DC'
vries2020@meta.data[vries2020@meta.data$cell_type %in% b, ]$cell_type_lowerres <- 'B'

# set cell type lowerres as the current identity
Idents(vries2020) <- 'cell_type_lowerres'

# read database of ligands and targets
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
# read ligand-receptor network
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
# read weighted integrated networks
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

# test for one receiver
expressed_genes_cd8t <- get_expressed_genes('CD8T', vries2020, pct = 0.10)
background_expressed_genes_cd8t <- expressed_genes_cd8t %>% .[. %in% rownames(ligand_target_matrix)]
# and multiple receivers
list_expressed_genes_notcd8t <- c('B', 'CD4T', 'DC', 'monocyte', 'NK') %>% unique() %>% lapply(get_expressed_genes, vries2020, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_notcd8t <- list_expressed_genes_notcd8t %>% unlist() %>% unique()

# extract CD8T
cd8t <- subset(vries2020, cell_type_lowerres == 'CD8T')
# set the condition as the ident
Idents(cd8t) <- 'stimulation'
# get the DE results from the receiving object
cd8t_de_table <- FindMarkers(cd8t, ident.1 = 'unstimulated', ident.2 = 'stimulated', min.pct = 0.1,  test.use = 'MAST')
# get the genes that were DE in cd8t stim/unstim
cd8t_de_table$gene <- rownames(cd8t_de_table)
cd8t_de_genes = cd8t_de_table %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
# filter agains ligand matrix
cd8t_de_genes = cd8t_de_genes %>% .[. %in% rownames(ligand_target_matrix)]

# get the unique ligands and receptors, pooling all cell types together
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

# get the expressed ligands and receptors
expressed_ligands_notcd8t = intersect(ligands,expressed_genes_notcd8t)
expressed_receptors_cd8t = intersect(receptors,expressed_genes_cd8t)

# do the actual analysis (need to look into more thoroughly)
ligand_activities = predict_ligand_activities(geneset = cd8t_de_genes, background_expressed_genes = background_expressed_genes_cd8t, ligand_target_matrix = ligand_target_matrix, potential_ligands = expressed_ligands_notcd8t)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

# get the top 20 genes
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

# 
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = cd8t_de_genes, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)
order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

# visualize
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))



