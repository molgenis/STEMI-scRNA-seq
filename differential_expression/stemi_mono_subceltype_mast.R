# Script for DE classical monocytes vs non-classical monocytes in scRNA-seq cardio samples

###################
# Libraries       #
###################

library(Seurat)
library(MAST)
library(UpSetR)

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

do_MAST_per_condition <- function(seurat_object, condition.1, condition.2, output_loc, condition.column='cell_type', split.column='timepoint.final', assay = 'RNA', min.pct = 0.1, features=NULL, logfc.threshold=0.25, latent.vars=NULL){
  # subset per timepoint combination
  seurat_object_ut <- seurat_object[, seurat_object[[split.column]] == 'UT']
  seurat_object_t24h <- seurat_object[, seurat_object[[split.column]] == 't24h' ]
  seurat_object_t8w <- seurat_object[, seurat_object[[split.column]] == 't8w' ]
  seurat_object_baseline <- seurat_object[, seurat_object[[split.column]] == 'Baseline']

  # call do_MAST with subsetted seurat object
  do_MAST(seurat_object_ut, paste(output_loc, 'UT_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_t24h, paste(output_loc, 't24h_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_t8w, output_loc, paste(output_loc, 't8w_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
  do_MAST(seurat_object_baseline, paste(output_loc, 'Baseline_', sep=''), condition.1 = condition.1, condition.2 = condition.2, assay = assay, min.pct = min.pct, features = features, logfc.threshold = logfc.threshold, latent.vars = latent.vars)
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

write_meta_mast <- function(mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend){
  # go through the conditions
  files <- list.files(paste(mast_output_loc_prepend, '2', mast_output_loc_append, sep = ''))
  for(file in files){
    # get the mast output
    mast_loc_v2 <- paste(mast_output_loc_prepend, '2', mast_output_loc_append, file, sep = '')
    # v3 should be the same, but then v3
    mast_loc_v3 <- paste(mast_output_loc_prepend, '3', mast_output_loc_append, gsub('v2', 'v3', file), sep = '')
    tryCatch({
      # read the mast output
      mast_v2 <- read.table(mast_loc_v2, header=T)
      mast_v3 <- read.table(mast_loc_v3, header=T)
      # get the genes that are in both
      genes_both <- intersect(rownames(mast_v2), rownames(mast_v3))
      # select only those genes
      mast_v2 <- mast_v2[rownames(mast_v2) %in% genes_both,]
      mast_v3 <- mast_v3[rownames(mast_v3) %in% genes_both,]
      # morph P val to minimum
      if(nrow(mast_v2[mast_v2$p_val == 0, ]) > 0){
        mast_v2[mast_v2$p_val == 0, ]$p_val <- .Machine$double.xmin
      }
      if(nrow(mast_v3[mast_v3$p_val == 0, ]) > 0){
        mast_v3[mast_v3$p_val == 0, ]$p_val <- .Machine$double.xmin
      }
      # add the gene name also in a column
      mast_v2$gene <- rownames(mast_v2)
      mast_v3$gene <- rownames(mast_v3)
      # add the mast results
      masts <- list()
      masts$v2 <- mast_v2
      masts$v3 <- mast_v3
      # perform the metavolcanor approach
      meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val', foldchangecol='avg.logFC', genenamecol = 'gene', collaps = T)
      # grab the result we care about
      volcanometa <- meta_degs_comb@metaresult
      volcanometa$metap_bonferroni <- volcanometa$metap*length(genes_both)
      # add the genes as rownames
      rownames(volcanometa) <- volcanometa$gene
      # add a colname append
      colnames(mast_v2) <- paste(colnames(mast_v2), 'v2', sep = '_')
      colnames(mast_v3) <- paste(colnames(mast_v3), 'v3', sep = '_')
      # merge the frames
      mast <- merge(mast_v2, mast_v3, by=0, all=TRUE)
      rownames(mast) <- mast$Row.names
      mast$Row.names <- NULL
      # also add the volcanometa stuff
      mast <- merge(mast, volcanometa, by=0, all=TRUE)
      rownames(mast) <- mast$Row.names
      mast$Row.names <- NULL
      if(nrow(mast[mast$metap_bonferroni > 1, ]) > 0){
        mast[mast$metap_bonferroni > 1, ]$metap_bonferroni <- 1
      }
      # write the result
      output_loc <- paste(mast_meta_output_loc_prepend, file, sep = '')
      write.table(mast, output_loc, sep = '\t')
    }, error = function(e) {
      print(paste('error in', file, 'due to', e))
    })
  }
}

get_sig_genes_mono_and_subceltypes <- function(condition_compare_loc, subceltype_compare_loc, condition.1, condition.2, pval_column='metap_bonferroni', sig_pval=0.05){
  # the file where the mono conditions are compared
  mono_conditions_output_loc <- paste(condition_compare_loc, 'monocyte', condition.1, condition.2, '.tsv', sep = '')
  print(mono_conditions_output_loc)
  # the file where mono1 and mono2 are compared for condition 1
  mono_subcell_condition1_output_loc <- paste(subceltype_compare_loc, condition.1, '_mono 1mono 2.tsv', sep = '')
  print(mono_subcell_condition1_output_loc)
  # the file where mono1 and mono2 are compared for condition 2
  mono_subcell_condition2_output_loc <- paste(subceltype_compare_loc, condition.2, '_mono 1mono 2.tsv', sep = '')
  print(mono_subcell_condition2_output_loc)
  # init the result
  result <- list()
  tryCatch({
    # read the mast output
    mast_conditions <- read.table(mono_conditions_output_loc, header=T)
    # filter to only include the significant results
    mast_conditions <- mast_conditions[mast_conditions[[pval_column]] <= sig_pval, ]
    # same for the separate conditions
    mast_condition1 <- read.table(mono_subcell_condition1_output_loc, header = T)
    mast_condition1 <- mast_condition1[mast_condition1[[pval_column]] <= sig_pval, ]
    mast_condition2 <- read.table(mono_subcell_condition2_output_loc, header = T)
    mast_condition2 <- mast_condition2[mast_condition2[[pval_column]] <= sig_pval, ]
    # put in in the list
    name_conditions <- paste(condition.1, condition.2, sep=' ')
    result[[name_conditions]] <- rownames(mast_conditions)
    name_condition1 <- paste(condition.1, 'c vs nc')
    result[[name_condition1]] <- rownames(mast_condition1)
    name_condition2 <- paste(condition.2, 'c vs nc')
    result[[name_condition2]] <- rownames(mast_condition2)
  }, error = function(e) {
    print(paste('error in', condition.1, ' and ', condition.2, 'due to', e))
  })
  return(result)
}


###################
# Main code       #
###################

mast_output_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/'
mast_output_loc_v2 <- paste(mast_output_loc, 'stemi_v2_cmono_ncmono_20200911/rna/', sep='')
mast_output_loc_v3 <- paste(mast_output_loc, 'stemi_v3_cmono_ncmono_20200911/rna/', sep='')

# load object
cardio.integrated = readRDS(file = "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.integrated_20200820.rds")
# do the work for v2
cardio.chem2 <- subset(cardio.integrated, subset = chem == 'V2')
DefaultAssay(cardio.chem2) <- 'RNA'
cardio.chem2 <- NormalizeData(cardio.chem2)
# call do_MAST_per_condition_combination with right parameters
do_MAST_per_condition_combination(cardio.chem2, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v2, logfc.threshold=0.1)
do_MAST_per_condition(cardio.chem2, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v2, logfc.threshold=0.1)
rm(cardio.chem2)

# do the work for v3
cardio.chem3 <- subset(cardio.integrated, subset = chem == 'V3')
DefaultAssay(cardio.chem3) <- 'RNA'
cardio.chem3 <- NormalizeData(cardio.chem3)
# call do_MAST_per_condition_combination with right parameters
do_MAST_per_condition_combination(cardio.chem3, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v3, logfc.threshold=0.1)
do_MAST_per_condition(cardio.chem3, condition.1 = 'mono 1', condition.2 = 'mono 2', output_loc = mast_output_loc_v3, logfc.threshold=0.1)
rm(cardio.chem3)


# set local locations now
mast_output_loc_prepend <- '/data/cardiology/differential_expression/MAST/results/stemi_v'
mast_output_loc_append <- '_cmono_ncmono_20200911/rna/'
mast_output_loc_meta <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_cmono_ncmono_20200911/rna/'

# load library here, as it's not working on the cluster
library(MetaVolcanoR)
write_meta_mast(mast_output_loc_prepend, mast_output_loc_append, mast_output_loc_meta)

# check cmono vs ncmono
upset(fromList(get_sig_genes_mono_and_subceltypes('/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc01minpct01ncountrna_20200707/rna/', mast_output_loc_meta, 'UT', 'Baseline')), order.by = 'freq')
upset(fromList(get_sig_genes_mono_and_subceltypes('/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc01minpct01ncountrna_20200707/rna/', mast_output_loc_meta, 'UT', 't24h')), order.by = 'freq')
upset(fromList(get_sig_genes_mono_and_subceltypes('/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc01minpct01ncountrna_20200707/rna/', mast_output_loc_meta, 'Baseline', 't24h')), order.by = 'freq')


