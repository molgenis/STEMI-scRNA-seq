# read the integrated object
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio_integrated_20pcs.rds')
# read stemi_v2
stemi_v2 <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_final_wdemuxcorrectedassignments.rds')
# read stemi_v3
stemi_v3 <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/stemi_v3_35pcs.rds')
# read the oneM eqtlgen object
oneM <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/1M_UT_eqtlgen.rds')

# perform RNA normalization on cardio.integrated, we might need it later
DefaultAssay(cardio.integrated) <- 'RNA'
cardio.integrated <- NormalizeData(cardio.integrated)
# harmonize the previous cell types
cardio.integrated@meta.data[(is.na(cardio.integrated@meta.data$cell_type)),]$cell_type <- 'unclassified'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'mDCs',]$cell_type <- 'mDC'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Megakaryocytes',]$cell_type <- 'megakaryocyte'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Memory_CD4T',]$cell_type <- 'memory CD4T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Memory_CD8T',]$cell_type <- 'memory CD8T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'mono2',]$cell_type <- 'mono 2'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Naive_CD8T',]$cell_type <- 'naive CD8T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'Naive_CD4T',]$cell_type <- 'naive CD4T'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'pDCs',]$cell_type <- 'pDC'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'PlasmaB',]$cell_type <- 'plasma B'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'cMonocytes',]$cell_type <- 'mono 1'
cardio.integrated@meta.data[cardio.integrated@meta.data$cell_type == 'ncMonocytes',]$cell_type <- 'mono 2'

# now store this in a new column to make sure we don't lose it
cardio.integrated@meta.data$cell_type_previous <- cardio.integrated@meta.data$cell_type
# fix this little inconsistency
levels(cardio.integrated@meta.data$orig.ident) <- c(levels(cardio.integrated@meta.data$orig.ident), 'stemi_v3')
cardio.integrated@meta.data[cardio.integrated@meta.data$orig.ident == 'cardio_v3',]$orig.ident <- 'stemi_v3'

# add the initial cell types
initial_ct <-  c(
  'mono 1',
  'naive CD4T',
  'memory CD8T',
  'NKdim', 
  'naive CD4T',
  'memory CD8T', 
  'memory CD4T',
  'mono 1',
  'memory CD4T',
  'memory CD4T',
  'naive CD4T',
  'mono 2',
  'NKdim',
  'B',
  'memory CD8T', 
  'mono 1',
  'memory CD4T',
  'mDC',
  'memory CD4T',
  'megakaryocyte',
  'pDC',
  'memory CD4T',
  'mono',
  'NKbright',
  'plasma B',
  'hemapoietic stem',
  'NKbright'
)
names(initial_ct) <- levels(cardio.integrated)
cardio.integrated <-RenameIdents(cardio.integrated, initial_ct)
# put it in a column as well
cardio.integrated <- AddMetaData(cardio.integrated, cardio.integrated@active.ident, 'cell_type')
cardio.integrated <- AddMetaData(cardio.integrated, cardio.integrated@active.ident, 'cell_type_initial')

# now save the monocyte celltyping
monocyte_split.integrated_15pcs = readRDS("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/monocyte_split.integrated_15pcs.rds")
# Normalising data
DefaultAssay(monocyte_split.integrated_15pcs) <- "RNA"
monocyte_split.integrated_15pcs <- NormalizeData(monocyte_split.integrated_15pcs)
# Naming/specifying clusters
monocyte_ct <- c("mono 1", "mono 1", "mono 1", "mono 1", "mono 1", "mono 2", "mono 1", "mono 1", "mono 1", "mono 3", "mono 1", "mono 1", "mono 1", "megakaryocyte", "mono 3", "mono 2", "mono 4", "mono 2", "mono 1", "mono 2")
names(monocyte_ct) <- levels(monocyte_split.integrated_15pcs)
monocyte_split.integrated_15pcs <- RenameIdents(monocyte_split.integrated_15pcs, monocyte_ct)
monocyte_split.integrated_15pcs <- AddMetaData(monocyte_split.integrated_15pcs, monocyte_split.integrated_15pcs@active.ident, "cell_type")
# already done -> saveRDS(monocyte_split.integrated_15pcs, "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/monocyte_integrated_15pcs_ctd_RNA.rds")
# grab the assignments from the object
mono_specific_assignments <- monocyte_split.integrated_15pcs@meta.data['cell_type']

# now save the T subtyping
t_split.integrated_15pcs <- readRDS("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/t_split.integrated_UMAP_15pcs.rds")
# Normalising data
DefaultAssay(t_split.integrated_15pcs) <- "RNA"
t_split.integrated_15pcs <- NormalizeData(t_split.integrated_15pcs)
# Naming/specifying clusters
#t_ct <- c('th1 memory CD4T', 'naive CD4T', 'naive CD4T', 'memory CD8T', 'naive CD8T', 'th2 memory CD4T', 'memory CD8T', 'memory CD4T', 'naive CD4T transitioning to stim', 'reg CD4T', 'th1 memory CD4T', 'memory CD8T', 'th2 memory CD4T', 'memory CD8T', 'Tc17 (IL17+ CD8+)', 'naive CD4T', 'naive CD4T', 'memory CD8T', 'megakaryocyte', 'naive CD4T', 'cyto memory CD4T')
t_ct <- c(
  'memory CD4T',
  'naive CD4T',
  'naive CD4T',
  'memory CD8T',
  'naive CD8T',
  'memory CD4T',
  'memory CD8T',
  'memory CD4T',
  'naive CD4T',
  'Treg', 
  'memory CD4T', 
  'memory CD8T', 
  'memory CD4T', 
  'memory CD8T',
  'memory CD8T',
  'naive CD4T',
  'naive CD4T',
  'memory CD8T', 
  'megakaryocyte',
  'naive CD4T', 
  'memory CD4T'
)
names(t_ct) <- levels(t_split.integrated_15pcs)
t_split.integrated_15pcs <- RenameIdents(t_split.integrated_15pcs, t_ct)
t_split.integrated_15pcs <- AddMetaData(t_split.integrated_15pcs, t_split.integrated_15pcs@active.ident, "cell_type")
# already done -> saveRDS(t_split.integrated_15pcs, "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/t_integrated_15pcs_ctd_RNA.rds")
# grab the assignments from the object
t_specific_assignments <- t_split.integrated_15pcs@meta.data['cell_type']


# add the T and mono
cardio.integrated <- AddMetaData(cardio.integrated, t_split.integrated_15pcs@meta.data['cell_type'], 'cell_type_t')
rm(t.integrated_15pcs)
cardio.integrated <- AddMetaData(cardio.integrated, monocyte_split.integrated_15pcs@meta.data['cell_type'], 'cell_type_mono')
rm(monocyte_split.integrated_15pcs)
levels(cardio.integrated@meta.data$cell_type) <- c(levels(cardio.integrated@meta.data$cell_type), levels(cardio.integrated@meta.data$cell_type_mono), levels(cardio.integrated@meta.data$cell_type_t))
cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type_t),]$cell_type <- cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type_t),]$cell_type_t
cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type_mono),]$cell_type <- cardio.integrated@meta.data[!is.na(cardio.integrated@meta.data$cell_type_mono),]$cell_type_mono
# save the object
saveRDS(cardio.integrated, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio_integrated_20pcs_ctTmono_20200620.rds')

