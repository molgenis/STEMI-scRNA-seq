### Merging selected Th17 cells in STEMI v2 and v3
#############################################################################
# Load libraries and object
#############################################################################
library(Seurat)
library(ggplot2)
library(MAST)
library(tidyr)

# Reading object
stemi_v2 = readRDS(file = "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_final_wdemuxcorrectedassignments.rds")
stemi_v3 = readRDS(file = "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v3_35pcs.rds")

# Looking for the names in object
colnames(stemi_v2@meta.data)
unique(stemi_v2@meta.data$chem)

colnames(stemi_v3@meta.data)
unique(stemi_v3@meta.data$chem)

# Using RNA default per chemistry
DefaultAssay(stemi_v2) <- "RNA"
stemi_v2 <- NormalizeData(stemi_v2)
DefaultAssay(stemi_v3) <- "RNA"
stemi_v3 <- NormalizeData(stemi_v3)

#############################################################################
# 10X Chemistry STEMI V2 Th17 and Tc17 subset definition ## Th17 was already defined by Monique van der Wijst (redoing this anyway, right?)
#############################################################################
foxplot_CCL20 <- FeaturePlot(stemi_v2, features = c("CCL20"))
foxplot_CD3E <- FeaturePlot(stemi_v2, features = c("CD3E"))

# Adding a metadata column for CCL20
stemi_v2 <- AddMetaData(stemi_v2, foxplot_CCL20$data["CCL20"], "CCL20_relative_expression")
stemi_v2@meta.data$above_CCL20_threshold <- F
stemi_v2@meta.data[stemi_v2@meta.data$CCL20_relative_expression > 1, ]$above_CCL20_threshold <- T

# Adding a metadata column for CD3E
stemi_v2 <- AddMetaData(stemi_v2, foxplot_CD3E$data["CD3E"], "CD3E_relative_expression")
stemi_v2@meta.data$above_CD3E_threshold <- F
stemi_v2@meta.data[stemi_v2@meta.data$CD3E_relative_expression > 1, ]$above_CD3E_threshold <- T

###### Combining cells ######
# Combining cells that are CCL20 and CD3E positive
stemi_v2@meta.data$above_CCL20_and_CD3E_threshold <- F
stemi_v2@meta.data[stemi_v2@meta.data$above_CCL20_threshold == T & stemi_v2@meta.data$above_CD3E_threshold == T , ]$above_CCL20_and_CD3E_threshold <- T

# Making a subset of the CCL20 and CD3E positive cells
stemi_v2_CCL20_CD3E <- subset(stemi_v2, subset = above_CCL20_and_CD3E_threshold == T)

# Selecting cluster 0 for Th17 and finding markers of selected cells vs. rest of cluster
stemi_v2_clus0 <- subset(stemi_v2, subset = seurat_clusters == 0)
Idents(stemi_v2_clus0) <- stemi_v2_clus0@meta.data$above_CCL20_and_CD3E_threshold
stemi_v2_clus0_thresholded_table <- FindMarkers(stemi_v2_clus0, ident.1=T, ident.2= F, test.use="MAST")
write.csv(stemi_v2_clus0_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v2_clus0_thresholded.csv")

# Findmarkers cluster 0 vs all other cells
Idents(stemi_v2) <- "seurat_clusters"
stemi_v2_clus0_thresholded_table <- FindMarkers(stemi_v2, ident.1=0, ident.2=NULL, test.use="MAST")
write.csv(stemi_v2_clus0_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v2_clus0_thresholded2.csv")

# Selecting cluster 9 for Tc17 and finding markers of selected cells vs. rest of cluster
stemi_v2_clus9 <- subset(stemi_v2, subset = seurat_clusters == 9)
Idents(stemi_v2_clus9) <- stemi_v2_clus9@meta.data$above_CCL20_and_CD3E_threshold
stemi_v2_clus0_thresholded_table <- FindMarkers(stemi_v2_clus9, ident.1=T, ident.2= F, test.use="MAST")
write.csv(stemi_v2_clus0_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v2_clus9_thresholded.csv")

# Findmarkers cluster 9 vs all other cells
Idents(stemi_v2) <- "seurat_clusters"
stemi_v2_clus9_thresholded_table <- FindMarkers(stemi_v2, ident.1=9, ident.2=NULL, test.use="MAST")
write.csv(stemi_v2_clus9_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v2_clus9_thresholded2.csv")

#############################################################################
# 10X Chemistry STEMI V3 Th17 and Tc17 subset definition
#############################################################################
foxplot_CCL20 <- FeaturePlot(stemi_v3, features = c("CCL20"))
foxplot_CD3E <- FeaturePlot(stemi_v3, features = c("CD3E"))

# Adding a metadata column for CCL20
stemi_v3 <- AddMetaData(stemi_v3, foxplot_CCL20$data["CCL20"], "CCL20_relative_expression")
stemi_v3@meta.data$above_CCL20_threshold <- F
stemi_v3@meta.data[stemi_v3@meta.data$CCL20_relative_expression > 1, ]$above_CCL20_threshold <- T

# Adding a metadata column for CD3E
stemi_v3 <- AddMetaData(stemi_v3, foxplot_CD3E$data["CD3E"], "CD3E_relative_expression")
stemi_v3@meta.data$above_CD3E_threshold <- F
stemi_v3@meta.data[stemi_v3@meta.data$CD3E_relative_expression > 1, ]$above_CD3E_threshold <- T

###### Combining cells ######
# Combining cells that are CCL20 and CD3E positive
stemi_v3@meta.data$above_CCL20_and_CD3E_threshold <- F
stemi_v3@meta.data[stemi_v3@meta.data$above_CCL20_threshold == T & stemi_v3@meta.data$above_CD3E_threshold == T , ]$above_CCL20_and_CD3E_threshold <- T

# Making a subset of the CCL20 and CD3E positive cells
stemi_v3_CCL20_CD3E <- subset(stemi_v3, subset = above_CCL20_and_CD3E_threshold == T)

# Selecting cluster 0 for Th17 and finding markers of selected cells vs. rest of cluster
stemi_v3_clus0 <- subset(stemi_v3, subset = seurat_clusters == 0)
Idents(stemi_v3_clus0) <- stemi_v3_clus0@meta.data$above_CCL20_and_CD3E_threshold
stemi_v3_clus0_thresholded_table <- FindMarkers(stemi_v3_clus0, ident.1=T, ident.2= F, test.use="MAST")
write.csv(stemi_v3_clus0_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_clus0_thresholded.csv")

# Findmarkers cluster 0 (Th17) vs all other cells
Idents(stemi_v3) <- "seurat_clusters"
stemi_v3_clus0_thresholded_table <- FindMarkers(stemi_v3, ident.1=0, ident.2=NULL, test.use="MAST")
write.csv(stemi_v3_clus0_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_clus0_thresholded2.csv")

# Selecting cluster 9 for Tc17 and finding markers of selected cells vs. rest of cluster
stemi_v3_clus9 <- subset(stemi_v3, subset = seurat_clusters == 9)
Idents(stemi_v3_clus9) <- stemi_v3_clus9@meta.data$above_CCL20_and_CD3E_threshold
stemi_v3_clus9_thresholded_table <- FindMarkers(stemi_v3_clus9, ident.1=T, ident.2= F, test.use="MAST")
write.csv(stemi_v3_clus9_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_clus9_thresholded.csv")

# Findmarkers cluster 9 (Tc17) vs all other cells
Idents(stemi_v3) <- "seurat_clusters"
stemi_v3_clus9_thresholded_table <- FindMarkers(stemi_v3, ident.1=9, ident.2=NULL, test.use="MAST")
write.csv(stemi_v3_clus9_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_clus9_thresholded2.csv")

# Findmarkers of cluster 0 and 9 (Finding IL4I1)
stemi_v3_clus09 <- subset(stemi_v3, subset = seurat_clusters == 0 | seurat_clusters == 9)
Idents(stemi_v3_clus09) <- stemi_v3_clus09@meta.data$above_CCL20_and_CD3E_threshold
stemi_v3_clus09_thresholded_table <- FindMarkers(stemi_v3_clus09, ident.1=T, ident.2= F, test.use="MAST")
write.csv(stemi_v3_clus09_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_clus09_thresholded.csv")

#############################################################################
# Adding STEMI V2 Th17 and Tc17 subset to V2 object
#############################################################################
# create new column
stemi_v2@meta.data$cell_type_better <- as.factor('unknown') #NOTE THAT THIS DOES NOT TAKE Th17 INTO ACCOUNT ALREADY IN THE CELL_TYPES COLUMN!
# add new options to column
levels(stemi_v2@meta.data$cell_type_better) <- c(levels(stemi_v2@meta.data$cell_type_better), "Th17", "Tc17")
# set new cell types
stemi_v2@meta.data[stemi_v2@meta.data$above_CCL20_and_CD3E_threshold ==T & stemi_v2@meta.data$seurat_clusters == 0 ,]$cell_type_better <- "Th17"
stemi_v2@meta.data[stemi_v2@meta.data$above_CCL20_and_CD3E_threshold ==T & stemi_v2@meta.data$seurat_clusters == 9 ,]$cell_type_better <- "Tc17"

#############################################################################
# Adding STEMI V3 Th17 and Tc17 subset to V3 object
#############################################################################
# create new column
stemi_v3@meta.data$cell_type_better <- as.factor('unknown')
# add new options to column
levels(stemi_v3@meta.data$cell_type_better) <- c(levels(stemi_v3@meta.data$cell_type_better), "Th17", "Tc17")
# set new cell types
stemi_v3@meta.data[stemi_v3@meta.data$above_CCL20_and_CD3E_threshold ==T & stemi_v3@meta.data$seurat_clusters == 0 ,]$cell_type_better <- "Th17"
stemi_v3@meta.data[stemi_v3@meta.data$above_CCL20_and_CD3E_threshold ==T & stemi_v3@meta.data$seurat_clusters == 9 ,]$cell_type_better <- "Tc17"

#############################################################################
# Plotting Th17 and Tc17 in STEMI V2 and V3
#############################################################################
# Plotting previously defined cell types
DimPlot(stemi_v2, group.by= "cell_type_better")
ggsave("/home/umcg-ivanblokland/stemi/stemi_v2_celltypeTh17.png", dpi=600, height=20, width=20)
DimPlot(stemi_v3, group.by= "cell_type_better")
ggsave("/home/umcg-ivanblokland/stemi/stemi_v3_celltypeTh17.png", dpi=600, height=20, width=20)

# Plotting Th17 and Tc17 <- beed to add initial cell types first!
#DimPlot(stemi_v2, group.by= "cell_type_better", cols=c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey", "pink","orange"))
#ggsave("/home/umcg-ivanblokland/stemi/stemi_v2_celltypeTh17coloured.png", dpi=600, height=20, width=20)
#DimPlot(stemi_v3, group.by= "cell_type_better", cols=c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey", "pink","orange"))
#ggsave("/home/umcg-ivanblokland/stemi/stemi_v3_celltypeTh17coloured.png", dpi=600, height=20, width=20)

#############################################################################
# Counting amount of cells Th17 and Tc17 in STEMI V2 and V3
#############################################################################
# Counting Th17 cells
nrow(stemi_v2@meta.data[stemi_v2@meta.data$cell_type_better == "Th17",])
# [1] 47
nrow(stemi_v3@meta.data[stemi_v3@meta.data$cell_type_better == "Th17",])
# [1] 43

# Counting Tc17 cells
nrow(stemi_v2@meta.data[stemi_v2@meta.data$cell_type_better == "Tc17",])
# [1] 70
nrow(stemi_v3@meta.data[stemi_v3@meta.data$cell_type_better == "Tc17",])
# [1] 19

#############################################################################
# Checking STEMI V2 and V3 DE Th17 genes vs Tc17
#############################################################################
# V2: Findmarkers of Th17 vs Tc17
Idents(stemi_v2) <- "cell_type_better"
stemi_v2_Th17_thresholded_table <- FindMarkers(stemi_v2, ident.1="Th17", ident.2="Tc17", test.use="MAST")
write.csv(stemi_v2_Th17_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v2_Tc17vsTh17.csv")

# V3: Findmarkers of Th17 vs Tc17
Idents(stemi_v3) <- "cell_type_better"
stemi_v3_Th17_thresholded_table <- FindMarkers(stemi_v3, ident.1="Th17", ident.2="Tc17", test.use="MAST")
write.csv(stemi_v3_Th17_thresholded_table, file = "/home/umcg-ivanblokland/stemi/stemi_v3_Tc17vsTh17.csv")

#############################################################################
# Load libraries and object Healthy Controls (Untreated, UT)
#############################################################################
# Reading object
oneM_UT_eqtlgen <- readRDS(file = "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/1M_UT_eqtlgen.rds")

# Looking for the names in object
colnames(oneM_UT_eqtlgen@meta.data)
unique(oneM_UT_eqtlgen@meta.data$chem)

# Making separate objects for ut_1m_v2 and ut_1m_v3
ut_1m_v2 <- subset(oneM_UT_eqtlgen, subset = chem == "V2")
ut_1m_v3 <- subset(oneM_UT_eqtlgen, subset = chem == "V3")

# Using RNA default per chemistry
DefaultAssay(ut_1m_v2) <- "RNA"
ut_1m_v2 <- NormalizeData(ut_1m_v2)
DefaultAssay(ut_1m_v3) <- "RNA"
ut_1m_v3 <- NormalizeData(ut_1m_v3)

#############################################################################
# Making plots of UT V2 and V3
#############################################################################
# Overview V2
Idents(ut_1m_v2) <- "seurat_clusters"
DimPlot(ut_1m_v2, group.by = "seurat_clusters", label=T)
ggsave("/home/umcg-ivanblokland/stemi/oneM_UT_eqtlgen_v2_overview.png", dpi=600, height=20, width=20)

# Overview V3
Idents(ut_1m_v3) <- "seurat_clusters"
DimPlot(ut_1m_v3, group.by = "seurat_clusters", label=T)
ggsave("/home/umcg-ivanblokland/stemi/oneM_UT_eqtlgen_v3_overview.png", dpi=600, height=20, width=20)

#############################################################################
# UT V2 Th17 subset definition
#############################################################################
foxplot_CCL20 <- FeaturePlot(ut_1m_v2, features = c("CCL20"))
foxplot_CD3E <- FeaturePlot(ut_1m_v2, features = c("CD3E"))

# Adding a metadata column for CCL20
ut_1m_v2 <- AddMetaData(ut_1m_v2, foxplot_CCL20$data["CCL20"], "CCL20_relative_expression")
ut_1m_v2@meta.data$above_CCL20_threshold <- F
ut_1m_v2@meta.data[ut_1m_v2@meta.data$CCL20_relative_expression > 1, ]$above_CCL20_threshold <- T

# Adding a metadata column for CD3E
ut_1m_v2 <- AddMetaData(ut_1m_v2, foxplot_CD3E$data["CD3E"], "CD3E_relative_expression")
ut_1m_v2@meta.data$above_CD3E_threshold <- F
ut_1m_v2@meta.data[ut_1m_v2@meta.data$CD3E_relative_expression > 1, ]$above_CD3E_threshold <- T

###### Combining cells ######
# Combining cells that are CCL20 and CD3E positive
ut_1m_v2@meta.data$above_CCL20_and_CD3E_threshold <- F
ut_1m_v2@meta.data[ut_1m_v2@meta.data$above_CCL20_threshold == T & ut_1m_v2@meta.data$above_CD3E_threshold == T , ]$above_CCL20_and_CD3E_threshold <- T

# Making a subset of the CCL20 and CD3E positive cells
ut_1m_v2_CCL20_CD3E <- subset(ut_1m_v2, subset = above_CCL20_and_CD3E_threshold == T)

# Selecting cluster 2 for Th17 and finding markers of selected cells vs. rest of cluster
ut_1m_v2_clus2 <- subset(ut_1m_v2, subset = seurat_clusters == 2)
Idents(ut_1m_v2_clus2) <- ut_1m_v2_clus2@meta.data$above_CCL20_and_CD3E_threshold
ut_1m_v2_clus2_thresholded_table <- FindMarkers(ut_1m_v2_clus2, ident.1=T, ident.2= F, test.use="MAST")
write.csv(ut_1m_v2_clus2_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v2_clus2_thresholded.csv")

# Findmarkers cluster 2 (Th17) vs all other cells
Idents(ut_1m_v2) <- "seurat_clusters"
ut_1m_v2_clus2_thresholded_table <- FindMarkers(ut_1m_v2, ident.1=2, ident.2=NULL, test.use="MAST")
write.csv(ut_1m_v2_clus2_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v2_clus2_thresholded2.csv")

# Selecting cluster 1 for Tc17 and finding markers of selected cells vs. rest of cluster
ut_1m_v2_clus1 <- subset(ut_1m_v2, subset = seurat_clusters == 1)
Idents(ut_1m_v2_clus1) <- ut_1m_v2_clus1@meta.data$above_CCL20_and_CD3E_threshold
ut_1m_v2_clus1_thresholded_table <- FindMarkers(ut_1m_v2_clus1, ident.1=T, ident.2= F, test.use="MAST")
write.csv(ut_1m_v2_clus1_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v2_clus1_thresholded.csv")

# Findmarkers cluster 1 vs all other cells
Idents(ut_1m_v2) <- "seurat_clusters"
ut_1m_v2_clus1_thresholded_table <- FindMarkers(ut_1m_v2, ident.1=1, ident.2=NULL, test.use="MAST")
write.csv(ut_1m_v2_clus1_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v2_clus1_thresholded2.csv")

#############################################################################
# UT V3 Th17 subset definition
#############################################################################
foxplot_CCL20 <- FeaturePlot(ut_1m_v3, features = c("CCL20"))
foxplot_CD3E <- FeaturePlot(ut_1m_v3, features = c("CD3E"))

# Adding a metadata column for CCL20
ut_1m_v3 <- AddMetaData(ut_1m_v3, foxplot_CCL20$data["CCL20"], "CCL20_relative_expression")
ut_1m_v3@meta.data$above_CCL20_threshold <- F
ut_1m_v3@meta.data[ut_1m_v3@meta.data$CCL20_relative_expression > 1, ]$above_CCL20_threshold <- T

# Adding a metadata column for CD3E
ut_1m_v3 <- AddMetaData(ut_1m_v3, foxplot_CD3E$data["CD3E"], "CD3E_relative_expression")
ut_1m_v3@meta.data$above_CD3E_threshold <- F
ut_1m_v3@meta.data[ut_1m_v3@meta.data$CD3E_relative_expression > 1, ]$above_CD3E_threshold <- T

###### Combining cells ######
# Combining cells that are CCL20 and CD3E positive
ut_1m_v3@meta.data$above_CCL20_and_CD3E_threshold <- F
ut_1m_v3@meta.data[ut_1m_v3@meta.data$above_CCL20_threshold == T & ut_1m_v3@meta.data$above_CD3E_threshold == T , ]$above_CCL20_and_CD3E_threshold <- T

# Making a subset of the CCL20 and CD3E positive cells
ut_1m_v3_CCL20_CD3E <- subset(ut_1m_v3, subset = above_CCL20_and_CD3E_threshold == T)

# Selecting cluster 2 as Th17 and finding markers of selected cells vs. rest of cluster
ut_1m_v3_clus2 <- subset(ut_1m_v3, subset = seurat_clusters == 2)
Idents(ut_1m_v3_clus2) <- ut_1m_v3_clus2@meta.data$above_CCL20_and_CD3E_threshold
ut_1m_v3_clus2_thresholded_table <- FindMarkers(ut_1m_v3_clus2, ident.1=T, ident.2= F, test.use="MAST")
write.csv(ut_1m_v3_clus2_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v3_clus2_thresholded.csv")

# Findmarkers cluster 2 (Th17) vs all other cells
Idents(ut_1m_v3) <- "seurat_clusters"
ut_1m_v3_clus2_thresholded_table <- FindMarkers(ut_1m_v3, ident.1=2, ident.2=NULL, test.use="MAST")
write.csv(ut_1m_v3_clus2_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v3_clus2_thresholded2.csv")

# Selecting cluster 1 for Tc17 and finding markers of selected cells vs. rest of cluster
ut_1m_v3_clus1 <- subset(ut_1m_v3, subset = seurat_clusters == 1)
Idents(ut_1m_v3_clus1) <- ut_1m_v3_clus1@meta.data$above_CCL20_and_CD3E_threshold
ut_1m_v3_clus1_thresholded_table <- FindMarkers(ut_1m_v3_clus1, ident.1=T, ident.2= F, test.use="MAST")
write.csv(ut_1m_v3_clus1_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v3_clus1_thresholded.csv")

# Findmarkers cluster 1 vs all other cells
Idents(ut_1m_v3) <- "seurat_clusters"
ut_1m_v3_clus1_thresholded_table <- FindMarkers(ut_1m_v3, ident.1=1, ident.2=NULL, test.use="MAST")
write.csv(ut_1m_v3_clus1_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v3_clus1_thresholded2.csv")

#############################################################################
# Adding UT V2 Th17 and Tc17 subset to UT V2 object
#############################################################################
# create new column
ut_1m_v2@meta.data$cell_type_better <- ut_1m_v2@meta.data$cell_type
# add new options to column
levels(ut_1m_v2@meta.data$cell_type_better) <- c(levels(ut_1m_v2@meta.data$cell_type_better), "Th17", "Tc17")
# set new cell types
ut_1m_v2@meta.data[ut_1m_v2@meta.data$above_CCL20_and_CD3E_threshold ==T & ut_1m_v2@meta.data$seurat_clusters == 2 ,]$cell_type_better <- "Th17"
ut_1m_v2@meta.data[ut_1m_v2@meta.data$above_CCL20_and_CD3E_threshold ==T & ut_1m_v2@meta.data$seurat_clusters == 1 ,]$cell_type_better <- "Tc17"

#############################################################################
# Adding UT V3 Th17 and Tc17 subset to UT V3 object
#############################################################################
# create new column
ut_1m_v3@meta.data$cell_type_better <- ut_1m_v3@meta.data$cell_type
# add new options to column
levels(ut_1m_v3@meta.data$cell_type_better) <- c(levels(ut_1m_v3@meta.data$cell_type_better), "Th17", "Tc17")
# set new cell types
ut_1m_v3@meta.data[ut_1m_v3@meta.data$above_CCL20_and_CD3E_threshold ==T & ut_1m_v3@meta.data$seurat_clusters == 2 ,]$cell_type_better <- "Th17"
ut_1m_v3@meta.data[ut_1m_v3@meta.data$above_CCL20_and_CD3E_threshold ==T & ut_1m_v3@meta.data$seurat_clusters == 1 ,]$cell_type_better <- "Tc17"

#############################################################################
# Plotting Th17 and Tc17 in UT V2 and V3
#############################################################################
# Plotting previously defined cell types
DimPlot(ut_1m_v2, group.by= "cell_type_better")
ggsave("/home/umcg-ivanblokland/stemi/ut_1m_v2_celltypeTh17.png", dpi=600, height=20, width=20)
DimPlot(ut_1m_v3, group.by= "cell_type_better")
ggsave("/home/umcg-ivanblokland/stemi/ut_1m_v3_celltypeTh17.png", dpi=600, height=20, width=20)

# Plotting Th17 and Tc17
DimPlot(ut_1m_v2, group.by= "cell_type_better", cols=c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey", "pink","orange"))
ggsave("/home/umcg-ivanblokland/stemi/ut_1m_v2_celltypeTh17coloured.png", dpi=600, height=20, width=20)
DimPlot(ut_1m_v3, group.by= "cell_type_better", cols=c("grey","grey","grey","grey","grey","grey","grey","grey","grey","grey","grey", "pink","orange"))
ggsave("/home/umcg-ivanblokland/stemi/ut_1m_v3_celltypeTh17coloured.png", dpi=600, height=20, width=20)

#############################################################################
# Counting amount of cells Th17 and Tc17 in UT V2 and V3
#############################################################################
# Counting Th17 cells
nrow(ut_1m_v2@meta.data[ut_1m_v2@meta.data$cell_type_better == "Th17",])
# [1] 439
nrow(ut_1m_v3@meta.data[ut_1m_v3@meta.data$cell_type_better == "Th17",])
# [1] 304

# Counting Tc17 cells
nrow(ut_1m_v2@meta.data[ut_1m_v2@meta.data$cell_type_better == "Tc17",])
# [1] 310
nrow(ut_1m_v3@meta.data[ut_1m_v3@meta.data$cell_type_better == "Tc17",])
# [1] 571

#############################################################################
# Checking DE Th17 markers vs Tc17
#############################################################################
# V2: Findmarkers cluster 2 vs all other cells
Idents(ut_1m_v2) <- "cell_type_better"
ut_1m_v2_Th17_thresholded_table <- FindMarkers(ut_1m_v2, ident.1="Th17", ident.2="Tc17", test.use="MAST")
write.csv(ut_1m_v2_Th17_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v2_Tc17vsTh17.csv")

# V3: Findmarkers cluster 2 vs all other cells
Idents(ut_1m_v3) <- "cell_type_better"
ut_1m_v3_Th17_thresholded_table <- FindMarkers(ut_1m_v3, ident.1="Th17", ident.2="Tc17", test.use="MAST")
write.csv(ut_1m_v3_Th17_thresholded_table, file = "/home/umcg-ivanblokland/stemi/ut_1m_v3_Tc17vsTh17.csv")

#############################################################################
# Adding 
#############################################################################
# add the separate metadatas back together
remerged_metadata_ut_ct <- rbind(ut_1m_v2@meta.data['cell_type_better'], ut_1m_v3@meta.data['cell_type_better'])
# add the column we care about to the merged object
oneM_UT_eqtlgen <- AddMetaData(oneM_UT_eqtlgen, remerged_metadata_ut_ct["cell_type_better"], "cell_type_betterest")
# checking results
DimPlot(oneM_UT_eqtlgen, group.by = "cell_type_betterest")
ggsave("/home/umcg-ivanblokland/stemi/oneM_UT_eqtlgen_th17tc17.png", width = 10, height = 10, dpi = 600)

#############################################################################
# Saving the work
#############################################################################
# Saving
saveRDS(oneM_UT_eqtlgen, "/home/umcg-ivanblokland/stemi/oneM_UT_eqtlgen_new.rds")

#############################################################################
# Merging STEMI and UT
#############################################################################
# read integrated object
cardio.integrated <- readRDS('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio_integrated_20pcs_ctTmono_20200620.rds')

# backup old celltype
cardio.integrated@meta.data$cell_type_no17 <- cardio.integrated@meta.data$cell_type
# add new celltypes as options
levels(cardio.integrated@meta.data$cell_type) <- c(levels(cardio.integrated@meta.data$cell_type), 'Th17', 'Tc17')
# add Th17 from oneM
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(oneM_UT_eqtlgen@meta.data[oneM_UT_eqtlgen@meta.data$cell_type_betterest == 'Th17',]),]$cell_type <- 'Th17'
# add Tc17 from oneM
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(oneM_UT_eqtlgen@meta.data[oneM_UT_eqtlgen@meta.data$cell_type_betterest == 'Tc17',]),]$cell_type <- 'Tc17'
# add Th17 from stemi_v3
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(stemi_v3@meta.data[stemi_v3@meta.data$cell_type_better == 'Th17',]),]$cell_type <- 'Th17'
# add Tc17 from stemi_v3
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(stemi_v3@meta.data[stemi_v3@meta.data$cell_type_better == 'Tc17',]),]$cell_type <- 'Tc17'
# add Th17 from stemi_v2
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(stemi_v2@meta.data[stemi_v2@meta.data$cell_type_better == 'Th17',]),]$cell_type <- 'Th17'
# add Tc17 from stemi_v2
cardio.integrated@meta.data[rownames(cardio.integrated@meta.data) %in% rownames(stemi_v2@meta.data[stemi_v2@meta.data$cell_type_better == 'Tc17',]),]$cell_type <- 'Tc17'

# set cell type as default idents
Idents(cardio.integrated) <- 'cell_type'
# create dimplot
DimPlot(cardio.integrated)
ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.20pcs.ctTmonoTc17Th17.20200620.png', dpi = 600, width = 20, height = 20)
DimPlot(cardio.integrated, cols = c('gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'blue', 'red'))
ggsave('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/plots/dimplots/cardio.integrated.20pcs.ctTmonoTc17Th17.20200620_17color.png', dpi = 600, width = 20, height = 20)

# Counting Tc17 cells
nrow(cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == "Tc17",])
# [1] 404
nrow(cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == "Th17",])
# [1] 322

# save the end result
saveRDS(cardio.integrated, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio_integrated_20pcs_ctTmonoTc17Th17_20200620.rds')
