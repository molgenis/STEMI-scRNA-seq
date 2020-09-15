# Script for DE classical monocytes vs non-classical monocytes in scRNA-seq cardio samples

###################
# Libraries       #
###################

library(Seurat)
library(MAST)

###################
# Functions       #
###################

do_MAST_per_condition_combination <- function(seurat_object, condition.1, condition.2, output_loc, condition.column='cell_type', split.column='timepoint.final', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.25, latent.vars=NULL){
  # subset per timepoint combination
  seurat_object_ut_baseline <- seurat_object[, seurat_object[[split.column]] == 'UT' | seurat_object[[split.column]] == 'Baseline' ]
  seurat_object_ut_t24h <- seurat_object[, seurat_object[[split.column]] == 'UT' | seurat_object[[split.column]] == 't24h' ]
  seurat_object_ut_t8w <- seurat_object[, seurat_object[[split.column]] == 'UT' | seurat_object[[split.column]] == 't8w' ]
  seurat_object_baseline_t24h <- seurat_object[, seurat_object[[split.column]] == 'Baseline' | seurat_object[[split.column]] == 't24h' ]
  seurat_object_baseline_t8w <- seurat_object[, seurat_object[[split.column]] == 'Baseline' | seurat_object[[split.column]] == 't8w' ]
  seurat_object_t24h_t8w <- seurat_object[, seurat_object[[split.column]] == 't24h' | seurat_object[[split.column]] == 't8w' ]
  
  # call do_MAST with subsetted seurat object
  do_MAST(seurat_object_ut_baseline, paste(output_loc, 'UTBaseline_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_ut_t24h, paste(output_loc, 'UTt24h_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_ut_t8w, output_loc, paste(output_loc, 'UTt8w_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_baseline_t24h, paste(output_loc, 'baselinet24h', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_baseline_t8w, paste(output_loc, 'baselinet8w', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_t24h_t8w, paste(output_loc, 't24ht8w', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
}

do_MAST <- function(seurat_object, condition.1, condition.2, output_loc, condition.column='cell_type', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.25, latent.vars=NULL){
  # set the active ident to the condition column
  Idents(seurat_object) <- condition.column
  
  # Perform MAST if possible
  result <- NULL
  # call MAST with parameters in function call
  try({
    result <- FindMarkers(object = seurat_object, ident.1 = condition.1, ident.2 = condition.2, test.use = 'MAST', min.pct = min.pct, assay = assay, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  })

  # create an output directory
  output_loc_final <- paste(output_loc, condition.1, condition.2, '.tsv', sep = '')
  # write result somewhere
  if(is.null(result)){
    print(paste('nothing to compare due to threshold for', output_loc_final))
  }
  else{
    write.table(result, output_loc_final, sep = '\t')
  }
}

###################
# Main code       #
###################

mast_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
mast_output_loc_v2 <- paste(mast_output_loc, 'stemi_v2_cmono_ncmono_20200911/rna/')
mast_output_loc_v3 <- paste(mast_output_loc, 'stemi_v3_cmono_ncmono_20200911/rna/')

# load object
cardio.integrated = readRDS(file = "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated_20200820.rds")
# do the work for v2
cardio.chem2 <- subset(cardio.integrated, subset = chem == 'V2')
DefaultAssay(cardio.chem2) <- 'RNA'
cardio.chem2 <- NormalizeData(cardio.chem2)
# call do_MAST_per_condition_combination with right parameters
do_MAST_per_condition_combination(cardio.chem2, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v2, logfc.threshold=0.1)
rm(cardio.chem2)

# do the work for v3
cardio.chem3 <- subset(cardio.integrated, subset = chem == 'V3')
DefaultAssay(cardio.chem3) <- 'RNA'
cardio.chem3 <- NormalizeData(cardio.chem3)
# call do_MAST_per_condition_combination with right parameters
do_MAST_per_condition_combination(cardio.chem3, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v3, logfc.threshold=0.1)
rm(cardio.chem3)





