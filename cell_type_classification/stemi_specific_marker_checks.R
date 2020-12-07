library(Seurat)
library(ggplot2)


# plot celltype markers
plot_celltype_markers_l2 <- function(seurat_object, celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for 
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$predicted.celltype.l2 
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


# create violin plots
plot_celltype_violins_l2 <- function(seurat_object, celltype_marker_genes=c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34"), assay = "RNA", slot="data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # plot per gene
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = "predicted.celltype.l2", assay=assay, slot=slot)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10)
      VlnPlot(seurat_object, features = c(gene), group.by = "predicted.celltype.l2", assay=assay, slot=slot, pt.size = 0)
      ggsave(paste(plot_dir,gene,"_violin_nodot.png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
  # we don't want to error because of missing genes
  genes_to_plot <- intersect(rownames(seurat_object), celltype_marker_genes)
  # plot per cluster
  for(ident in levels(seurat_object@meta.data$predicted.celltype.l2)){
    # plot the violins and save
    VlnPlot(seurat_object, features = genes_to_plot, idents = c(ident), assay=assay, slot=slot)
    ggsave(paste(plot_dir,ident,"_violin.png", sep=""), dpi=600, width=25, height=25)
  }
}


# run findmarkers
run_findmarkers <- function(seurat_object, output_loc_final, ident_to_set, ident.1, ident.2=NULL){
  # set the assay we wish to analyse
  DefaultAssay(seurat_object) <- assay
  # a hard column name is required for the model.matrix step, so let's add that
  seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data[ident_to_set], 'condition')
  Idents(seurat_object) <- seurat_object@meta.data$condition
  # create an output directory
  #output_loc_final <- paste(output_loc, ident.1, ident.2, '.tsv', sep = '')
  # perform MAST if possible
  result <- NULL
  # we need to 'try' here, as MAST ungraciously exits if it cannot perform
  try({
    result <- FindMarkers(object = seurat_object, ident.1 = ident.1, ident.2 = ident.2, test.use = 'MAST', min.pct = min.pct, assay = assay, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  })
  # write the actual table if possible
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  } else{
    write.table(result, output_loc_final, sep = '\t')
  }
}

# set locations
# plots dir
plot_loc <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/"
# more specific plot locs
features_plot_loc <- paste(plot_loc, "feature_plots/", sep = "")
violins_plot_loc <- paste(plot_loc, "violin_plots/", sep = "")

# load the object
cardio.integrated <- readRDS('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated.20201126_wazi.rds')
# just to be sure, perform Normalization
DefaultAssay(cardio.integrated) <- 'RNA'
cardio.integrated <- NormalizeData(cardio.integrated)

# plot the cell type markers requested
plot_celltype_markers_l2(seurat_object = cardio.integrated, celltype_marker_genes=c('MKI67', 'HBB', 'HBA', 'HBC'), plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "irene/", sep = ""))
plot_celltype_violins_l2(seurat_object = cardio.integrated, celltype_marker_genes=c('MKI67', 'HBB', 'HBA', 'HBC'), plot_dir = paste(features_plot_loc, "cardio_integrated_20201126_30pcs/", "irene/", sep = ""))

# expression pattern is different in v2 and v3, better split these two up
cardio.integrated.v2 <- cardio.integrated[, cardio.integrated@meta.data$chem == 'V2']
cardio.integrated.v3 <- cardio.integrated[, cardio.integrated@meta.data$chem == 'V3']
# just to be sure, perform normalization
cardio.integrated.v2 <- NormalizeData(cardio.integrated.v2)
cardio.integrated.v3 <- NormalizeData(cardio.integrated.v3)
# clear memory
rm(cardio.integrated)

# subset to the CD4s
any_CD4 <- c('CD4 CTL', 'CD4 naive', 'CD4 TCM', 'CD4 TEM') # FILL THIS TO INCLUDE ALL possible CD4s!
cardio.integrated.v2.cd4 <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_CD4]
cardio.integrated.v3.cd4 <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_CD4]
# set the output location of our analyses regarding CD4
v2_cd4_ctl_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD4CTLvsCD4.tsv'
v2_cd4_tem_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD4TEMvsCD4.tsv'
v3_cd4_ctl_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD4CTLvsCD4.tsv'
v3_cd4_tem_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD4TEMvsCD4.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.cd4, output_loc_final=v2_cd4_ctl_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD4 CTL', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v2.cd4, output_loc_final=v2_cd4_tem_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD4 TEM', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd4, output_loc_final=v3_cd4_ctl_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD4 CTL', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd4, output_loc_final=v3_cd4_tem_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD4 TEM', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.cd4)
rm(cardio.integrated.v3.cd4)


# subset to CD8s
any_CD8 <- c('CD8 Naive', 'CD8 TCM', 'CD8 TEM')
cardio.integrated.v2.cd8 <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_CD8]
cardio.integrated.v3.cd8 <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_CD8]
# set the output location of our analyses regarding CD8
v2_cd8_naive_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD8NaivevsCD8.tsv'
v2_cd8_tcm_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD8TCMvsCD8.tsv'
v3_cd8_naive_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD8NaivevsCD8.tsv'
v3_cd8_tcm_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD8TCMvsCD8.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.cd8, output_loc_final=v2_cd8_naive_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD8 Naive', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v2.cd8, output_loc_final=v2_cd8_tcm_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD8 TCM', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd8, output_loc_final=v3_cd8_naive_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD8 Naive', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd8, output_loc_final=v3_cd8_tcm_out_loc, ident_to_set='predicted.celltype.l2', ident.1='CD8 TCM', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.cd8)
rm(cardio.integrated.v3.cd8)

# subset to CD8s
any_CD8_mait <- c('CD8 Naive', 'CD8 TCM', 'CD8 TEM', 'MAIT')
cardio.integrated.v2.cd8.mait <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_CD8_mait]
cardio.integrated.v3.cd8.mait <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_CD8_mait]
# set the output location of our analyses regarding CD8 MAIT
v2_cd8_mait_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD8MAITvsCD8.tsv'
v3_cd8_mait_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD8MAITvsCD8.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.cd8.mait, output_loc_final=v2_cd8_mait_out_loc, ident_to_set='predicted.celltype.l2', ident.1='MAIT', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd8.mait, output_loc_final=v3_cd8_mait_out_loc, ident_to_set='predicted.celltype.l2', ident.1='MAIT', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.cd8.mait)
rm(cardio.integrated.v3.cd8.mait)


# subset to CD8s
any_CD8_gdt <- c('CD8 Naive', 'CD8 TCM', 'CD8 TEM', 'gdT')
cardio.integrated.v2.cd8.gdt <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_CD8_gdt]
cardio.integrated.v3.cd8.gdt <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_CD8_gdt]
# set the output locaiton of our analysis regarding CD8 gdT
v2_cd8_gdt_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_CD8gdtvsCD8.tsv'
v3_cd8_gdt_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_CD8gdtvsCD8.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.cd8.gdt, output_loc_final=v2_cd8_gdt_out_loc, ident_to_set='predicted.celltype.l2', ident.1='gdT', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.cd8.gdt, output_loc_final=v3_cd8_gdt_out_loc, ident_to_set='predicted.celltype.l2', ident.1='gdT', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.cd8.gdt)
rm(cardio.integrated.v3.cd8.gdt)

# subset to NKs
#any_NK <- c('')
#cardio.integrated.v2.nk <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_NK]
#cardio.integrated.v3.nk <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_NK]
# set the output location of our analyses regarding NK
#v2_nk_prof_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_NKProvsNK.tsv'
#v3_nk_prof_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_NKProvsNK.tsv'

# subset to B
any_B <- c('B memory', 'B naive', 'B intermediate')
cardio.integrated.v2.b <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_B]
cardio.integrated.v3.b <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_B]
# set the output locatino of our analyses regarding B
v2_b_mem_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_BmemvsB.tsv'
v3_b_mem_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_BmemvsB.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.b, output_loc_final=v2_b_mem_out_loc, ident_to_set='predicted.celltype.l2', ident.1='B memory', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.b, output_loc_final=v3_b_mem_out_loc, ident_to_set='predicted.celltype.l2', ident.1='B memory', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.b)
rm(cardio.integrated.v3.b)

# subset to B with plasmablast
any_B_pb <- c('B memory', 'B naive', 'B intermediate', 'plasmablast')
cardio.integrated.v2.b.pb <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$predicted.celltype.l2 %in% any_B_pb]
cardio.integrated.v3.b.pb <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$predicted.celltype.l2 %in% any_B_pb]
# set the output locaiton of our analyses regarding B and plasmablasts
v2_b_plasblas_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_BplasblasvsB.tsv'
v3_b_plasblas_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_BplasblasvsB.tsv'
# do the actual analysis
run_findmarkers(seurat_object=cardio.integrated.v2.b, output_loc_final=v2_b_plasblas_out_loc, ident_to_set='predicted.celltype.l2', ident.1='plasmablast', ident.2=NULL)
run_findmarkers(seurat_object=cardio.integrated.v3.b, output_loc_final=v3_b_plasblas_out_loc, ident_to_set='predicted.celltype.l2', ident.1='plasmablast', ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.b.pb)
rm(cardio.integrated.v3.b.pb)

# subset to 1,4,8,13,17,25,31
cardio.integrated.v2.1.4.8.13.17.25.31 <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$seurat_clusters %in% c(1,4,8,13,17,25,31)]
cardio.integrated.v3.1.4.8.13.17.25.31 <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$seurat_clusters %in% c(1,4,8,13,17,25,31)]
# set the output location or analyses regarding 1,4,8,13,17,25,31
v2_1.4.8.13.17.25.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_25vs1.4.8.13.17.25.31.tsv'
v3_1.4.8.13.17.25.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_25vs1.4.8.13.17.25.31.tsv'
# do the actual analysis
run_findmarkers(cardio.integrated.v2.1.4.8.13.17.25.31, output_loc_final=v2_1.4.8.13.17.25.31_out_loc, ident_to_set='seurat_clusters', ident.1=25, ident.2=NULL)
run_findmarkers(cardio.integrated.v3.1.4.8.13.17.25.31, output_loc_final=v3_1.4.8.13.17.25.31_out_loc, ident_to_set='seurat_clusters', ident.1=25, ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.1.4.8.13.17.25.31)
rm(cardio.integrated.v3.1.4.8.13.17.25.31)


# subset to 1,4,8,13,31
cardio.integrated.v2.1.4.8.13.31 <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$seurat_clusters %in% c(1,4,8,13,31)]
cardio.integrated.v3.1.4.8.13.31 <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$seurat_clusters %in% c(1,4,8,13,31)]
# set the output location or of analyses regarding 1,4,8,13,31
v2_1.4.8.13.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_8vs1.4.8.13.31.tsv'
v3_1.4.8.13.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_8vs1.4.8.13.31.tsv'
# do the actual analysis
run_findmarkers(cardio.integrated.v2.1.4.8.13.31, output_loc_final=v2_1.4.8.13.31_out_loc, ident_to_set='seurat_clusters', ident.1=8, ident.2=NULL)
run_findmarkers(cardio.integrated.v3.1.4.8.13.31, output_loc_final=v3_1.4.8.13.31_out_loc, ident_to_set='seurat_clusters', ident.1=8, ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.1.4.8.13.31)
rm(cardio.integrated.v3.1.4.8.13.31)

# subset to 1,4,8,13,17,31
cardio.integrated.v2.1.4.8.13.17.31 <- cardio.integrated.v2[, cardio.integrated.v2@meta.data$seurat_clusters %in% c(1,4,8,13,17,31)]
cardio.integrated.v3.1.4.8.13.17.31 <- cardio.integrated.v3[, cardio.integrated.v3@meta.data$seurat_clusters %in% c(1,4,8,13,17,31)]
# set the output location or of analyses regarding 1,4,8,13,31
v2_1.4.8.13.17.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v2_17vs1.4.8.13.17.31.tsv'
v3_1.4.8.13.17.31_out_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/marker_output/cardio_integrated_20201126_30pcs/v3_17vs1.4.8.13.17.31.tsv'
# do the actual analysis
run_findmarkers(cardio.integrated.v2.1.4.8.13.17.31, output_loc_final=v2_1.4.8.13.17.31_out_loc, ident_to_set='seurat_clusters', ident.1=17, ident.2=NULL)
run_findmarkers(cardio.integrated.v3.1.4.8.13.17.31, output_loc_final=v3_1.4.8.13.17.31_out_loc, ident_to_set='seurat_clusters', ident.1=17, ident.2=NULL)
# clear memory
rm(cardio.integrated.v2.1.4.8.13.17.31)
rm(cardio.integrated.v3.1.4.8.13.17.31)
