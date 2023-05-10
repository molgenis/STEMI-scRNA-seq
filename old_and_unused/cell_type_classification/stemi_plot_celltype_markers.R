######################
# libraries          #
######################

library(Seurat)
library(ggplot2)

####################
# Functions        #
####################

# plot celltype markers
plot_celltype_markers <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters 
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# plot celltype markers
plot_celltype_markers_villani <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_villani <- c("CD1C","FCER1A","CLEC10A","CD1D","FCGR2B","CD33","S100A9","VCAN","LYZ","RNASE2","CD163","CSF3R","CD14","NAIP","F13A1","S100A12","CD36","TREM1","FCGR3A","LST1","AIF1","CTSS","MTSS1","TCF7L2","IFITM3","MS4A7","LILRB2","CSF1R","IFITM2","HCK","CTSL1","MAFB","TNFRSF1B","SIGLEC10","FGR","LILRA2","NEAT1","RHOC","EMR2","LAIR2","LILRB1","G0S2","NAMPT","FCGR3B","SRGN","TNFRSF10C","MXD1","CXCR2","VNN2","CXCR1","FCGR2A","CLEC4E","LITAF","TLR2","ITGB2","ITGAM","CTSD","CTSA","NLRP3","CLEC7A","BST1","STAB1","IRAK3","PRF1","GNLY","CTSW","KLRD1","NKG7","IL2RB","GZMA","ZAP70","KLRF1","GZMH","IL32","IKZF3","LCK","CD96","TGFBR3","CD2","CD247")
  celltype_marker_genes <- setdiff(celltype_marker_genes_villani, celltype_marker_genes_regular)
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# plot celltype markers
plot_celltype_markers_gate <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c("KRT1","NELL2","LMO7","ANKRD55","SLC37A3","SPINT2","EOMES","ARRDC5","RP11-326C3.11","PKIG","CCR6","CXCR6","B3GALT2","CCR5","PLEKHG3","LZTFL1","CCR2","KLRB1","CTLA4","GZMK","PTTG1","HMGN2","CCNB2","CAPG","KIAA0101","C12orf75","RANBP1","HNRNPA2B1","ANXA2","LMNA","HIST1H2AH","HIST1H2AJ","HIST1H1B","HIST1H2AM","TK1","CENPM","GTSE1","RRM2","HIST2H2AC","HIST1H1E","CCNB1","CDC20","HMMR","ASPM","CENPE","DLGAP4","UBE2C","CDKN3","KPNA2","CCL1","CCL3","MIR155HG","CCL4","FABP5","GZMB","TNFRSF4","IFNG","IL2RA","LTA","IL17RB","HPGDS","IL5","EGLN3","CREM","GADD45G","TNFSF10","GATA3","RPS10","PFDN5","FXYD5","LIMD2","SSR4","S100A11","PTPRCAP","CD3D","HLA-E","EMP3","FCER1G","RTKN2","LGALS3","RGS1","CCDC141","PIM1","PTGIR","CORO1B","TIGIT","C1orf56","APOBEC3C","NCR3","HNRNPH1","RASSF3_91","PPP1CB","SET","CDC42SE1","MDM4","B4GALT1")
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# plot celltype markers
plot_celltype_markers_gamez <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c("PLAC8","PTGER2","IL7R","NUC2B","HIST1H4C","HIST1H1B","HIST1H1E","TOP2A","GNRH2","EPAS1","MANF","TIGIT","GNLY","GZMA","CCL3","CCL4","NKG7","TMIGD2","STAG3","PRG4","PALLD","TMIGDIFIT3","IFIT1","RARRES3","RGS1","IL9","HLA-DRA","HLA-DRB1","BATF","CCL17","CD70","TNFRSF4","GZMH","BHLHE40","ANKRD37","HILPDA","DUSP4","IL2","DUSP2","REL","GATA3","MRPS26","LIMA1","SPINT2","CTLA4","TNFRSF18","IL2RA","BATF","IFIT2","PMAIP1","OASL","HERC5","DNAJB1","HSPA1A","HSPA1B","HSPB1","CCNB1","CDC20","CCNB2","UBE2S")
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# plot celltype markers
plot_celltype_markers_monique <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c('TBX21','ANXA3','GZMK','GATA3','MAOA','LIMA1','AUH','RUNX2','SMAD2','PALLD','BASP1','TNFRSF8','APOE','PIK3C2B','RUNX1','FOXP3','LGALS3','CCL5','GZMA','GZMB','IFIT1','IFI44L','IFI44','IFIT3','IFIT5','OAS3','OAS2','MX1','MX2','PRF1','GBP1','GBP4','SMAD9','CD96','CD38','IFNG','ALDOCS','AHR','CCL5','BASP1','RUNX1','ITGA4','CXCR4','BACH2','RRAS','FOXP4','PTK2','PRG4','CD70','TDRD7','TRIP10','RASGRP4','SEMA7A','LMNA','NCKAP1','ETHE1','LSR','TNFSF14','GZMB','CMPK2','CD38','STAT4','HERC6','ITGB7')
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# plot celltype markers
plot_celltype_markers_irene <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_irene <- c('FOXP3','CCR2','RORA','RORC','IL6','IL1','IL14','IL13','IL17A','SGP130','TREM2','VCAN ','CD16','CD16','IL-12','IL-23','IL-27','CXCL9','CXCL10','CXCL11','MMP1','TNFa','IL-10','CCL17')
  celltype_marker_genes <- setdiff(celltype_marker_genes_irene, celltype_marker_genes_regular)
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$seurat_clusters
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


# create violin plots
plot_celltype_violins <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}

# create violin plots
plot_celltype_violins_villani <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_villani <- c("CD1C","FCER1A","CLEC10A","CD1D","FCGR2B","CD33","S100A9","VCAN","LYZ","RNASE2","CD163","CSF3R","CD14","NAIP","F13A1","S100A12","CD36","TREM1","FCGR3A","LST1","AIF1","CTSS","MTSS1","TCF7L2","IFITM3","MS4A7","LILRB2","CSF1R","IFITM2","HCK","CTSL1","MAFB","TNFRSF1B","SIGLEC10","FGR","LILRA2","NEAT1","RHOC","EMR2","LAIR2","LILRB1","G0S2","NAMPT","FCGR3B","SRGN","TNFRSF10C","MXD1","CXCR2","VNN2","CXCR1","FCGR2A","CLEC4E","LITAF","TLR2","ITGB2","ITGAM","CTSD","CTSA","NLRP3","CLEC7A","BST1","STAB1","IRAK3","PRF1","GNLY","CTSW","KLRD1","NKG7","IL2RB","GZMA","ZAP70","KLRF1","GZMH","IL32","IKZF3","LCK","CD96","TGFBR3","CD2","CD247")
  celltype_marker_genes <- setdiff(celltype_marker_genes_villani, celltype_marker_genes_regular)
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}

# create violin plots
plot_celltype_violins_gate <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c("KRT1","NELL2","LMO7","ANKRD55","SLC37A3","SPINT2","EOMES","ARRDC5","RP11-326C3.11","PKIG","CCR6","CXCR6","B3GALT2","CCR5","PLEKHG3","LZTFL1","CCR2","KLRB1","CTLA4","GZMK","PTTG1","HMGN2","CCNB2","CAPG","KIAA0101","C12orf75","RANBP1","HNRNPA2B1","ANXA2","LMNA","HIST1H2AH","HIST1H2AJ","HIST1H1B","HIST1H2AM","TK1","CENPM","GTSE1","RRM2","HIST2H2AC","HIST1H1E","CCNB1","CDC20","HMMR","ASPM","CENPE","DLGAP4","UBE2C","CDKN3","KPNA2","CCL1","CCL3","MIR155HG","CCL4","FABP5","GZMB","TNFRSF4","IFNG","IL2RA","LTA","IL17RB","HPGDS","IL5","EGLN3","CREM","GADD45G","TNFSF10","GATA3","RPS10","PFDN5","FXYD5","LIMD2","SSR4","S100A11","PTPRCAP","CD3D","HLA-E","EMP3","FCER1G","RTKN2","LGALS3","RGS1","CCDC141","PIM1","PTGIR","CORO1B","TIGIT","C1orf56","APOBEC3C","NCR3","HNRNPH1","RASSF3_91","PPP1CB","SET","CDC42SE1","MDM4","B4GALT1")
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}

# create violin plots
plot_celltype_violins_gamez <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c("PLAC8","PTGER2","IL7R","NUC2B","HIST1H4C","HIST1H1B","HIST1H1E","TOP2A","GNRH2","EPAS1","MANF","TIGIT","GNLY","GZMA","CCL3","CCL4","NKG7","TMIGD2","STAG3","PRG4","PALLD","TMIGDIFIT3","IFIT1","RARRES3","RGS1","IL9","HLA-DRA","HLA-DRB1","BATF","CCL17","CD70","TNFRSF4","GZMH","BHLHE40","ANKRD37","HILPDA","DUSP4","IL2","DUSP2","REL","GATA3","MRPS26","LIMA1","SPINT2","CTLA4","TNFRSF18","IL2RA","BATF","IFIT2","PMAIP1","OASL","HERC5","DNAJB1","HSPA1A","HSPA1B","HSPB1","CCNB1","CDC20","CCNB2","UBE2S")
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}

# create violin plots
plot_celltype_violins_monique <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_gate <- c('TBX21','ANXA3','GZMK','GATA3','MAOA','LIMA1','AUH','RUNX2','SMAD2','PALLD','BASP1','TNFRSF8','APOE','PIK3C2B','RUNX1','FOXP3','LGALS3','CCL5','GZMA','GZMB','IFIT1','IFI44L','IFI44','IFIT3','IFIT5','OAS3','OAS2','MX1','MX2','PRF1','GBP1','GBP4','SMAD9','CD96','CD38','IFNG','ALDOCS','AHR','CCL5','BASP1','RUNX1','ITGA4','CXCR4','BACH2','RRAS','FOXP4','PTK2','PRG4','CD70','TDRD7','TRIP10','RASGRP4','SEMA7A','LMNA','NCKAP1','ETHE1','LSR','TNFSF14','GZMB','CMPK2','CD38','STAT4','HERC6','ITGB7')
  celltype_marker_genes <- setdiff(celltype_marker_genes_gate, celltype_marker_genes_regular)
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}

# create violin plots
plot_celltype_violins_irene <- function(seurat_object, assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes_regular <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  celltype_marker_genes_irene <- c('FOXP3','CCR2','RORA','RORC','IL6','IL1','IL14','IL13','IL17A','SGP130','TREM2','VCAN ','CD16','CD16','IL-12','IL-23','IL-27','CXCL9','CXCL10','CXCL11','MMP1','TNFa','IL-10','CCL17')
  celltype_marker_genes <- setdiff(celltype_marker_genes_irene, celltype_marker_genes_regular)
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "seurat_clusters", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$seurat_clusters)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}


#########################
# main code             #
#########################

# this is where our objects are on disk
object_loc <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
# objects specifically
cardio_30pcs_loc <- paste(object_loc, "cardio.integrated.20201126_wazi.rds", sep = "")

# plots dir
plot_loc <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/"
# more specific plot locs
features_plot_loc <- paste(plot_loc, "feature_plots/", sep = "")
violins_plot_loc <- paste(plot_loc, "violin_plots/", sep = "")

# do work for each object
cardio_30pcs <- readRDS(cardio_30pcs_loc)
# set to RNA and normalize (usual workflow for feature plots)
DefaultAssay(cardio_30pcs) <- "RNA"
cardio_30pcs <- NormalizeData(cardio_30pcs)
# plot various things
plot_celltype_markers(seurat_object = cardio_30pcs, plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "regular/", sep = ""))
plot_celltype_markers_gate(seurat_object = cardio_30pcs, plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "gate/", sep = ""))
plot_celltype_markers_villani(seurat_object = cardio_30pcs, plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "villani/", sep = ""))
plot_celltype_markers_irene(seurat_object = cardio_30pcs, plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "irene/", sep = ""))

plot_celltype_violins(seurat_object = cardio_30pcs, plot_dir = paste(violins_plot_loc, "cardio_integrated_20201126_30pcs/", "regular/", sep = ""))
plot_celltype_violins_gate(seurat_object = cardio_30pcs, plot_dir = paste(violins_plot_loc, "cardio_integrated_20201126_30pcs/", "gate/", sep = ""))
plot_celltype_violins_villani(seurat_object = cardio_30pcs, plot_dir = paste(violins_plot_loc, "cardio_integrated_20201126_30pcs/", "villani/", sep = ""))
plot_celltype_violins_irene(seurat_object = cardio_30pcs, plot_dir = paste(violins_plot_loc, "cardio_integrated_20201126_30pcs/", "irene/", sep = ""))

rm(cardio_30pcs)

