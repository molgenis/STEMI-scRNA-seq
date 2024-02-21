#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_anonimise_and_export_seurat.R
# Function: export the Seurat object into cellranger-like outputs
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(Matrix)

####################
# Functions        #
####################

create_anonymized_mapping <- function(metadata, source_column, target_column, target_prepend='') {
  # get the unique values
  source_values <- unique(as.character(metadata[[source_column]]))
  # turn NA into 'unknown'
  source_values[is.na(source_values)] <- 'unknown'
  # randomly sort your values, by sampling all values
  source_values <- source_values[sample(1:length(source_values), length(source_values))]
  # create a mapping
  target_values <- paste(target_prepend, 1:length(source_values), sep = '')
  # now put that into a dataframe
  mapping_table <- data.frame(x = source_values, y = target_values)
  # now make the column names as expected
  colnames(mapping_table) <- c(source_column, target_column)
  return(mapping_table)
}


####################
# Main Code        #
####################

# location of the seurat object
seurat_object_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds'
# location of the export
counts_export_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/export/'

# read the object
stemi <- readRDS(seurat_object_loc)

# features
features_gz <- gzfile(paste(counts_export_loc, 'features.tsv.gz', sep = ''))
write.table(data.frame(x = rownames(stemi@assays$SCT@counts)), features_gz, row.names = F, col.names = F, quote = F)
# barcodes
barcodes_gz <- gzfile(paste(counts_export_loc, 'barcodes.tsv.gz', sep = ''))
write.table(data.frame(x = colnames(stemi@assays$SCT@counts)), barcodes_gz, row.names = F, col.names = F, quote = F)
# and finally the count matrix
writeMM(stemi@assays$SCT@counts, paste(counts_export_loc, 'matrix.mtx', sep = ''))
# matrix was manually zipped with bgzip /groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/export/matrix.mtx

# unnormalized as well
features_gz_raw <- gzfile(paste(counts_export_loc, 'features_raw.tsv.gz', sep = ''))
write.table(data.frame(x = rownames(stemi@assays$RNA@counts)), features_gz_raw, row.names = F, col.names = F, quote = F)
writeMM(stemi@assays$RNA@counts, paste(counts_export_loc, 'matrix_raw.mtx', sep = ''))
# matrix was manually zipped with bgzip /groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/export/matrix.mtx

# get the sample names for the controls, as they were in the Oelen et. al. 2022 study
ll_mapping_loc <- '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/1M_exp_to_ll.tsv'
ll_mappping <- read.table(ll_mapping_loc, header = T)
# rename columns
colnames(ll_mappping) <- c('sample', 'geno_id')
# also add the 'C' to the sample
ll_mappping[['sample']] <- paste('C', ll_mappping[['sample']], sep = '')
# get a mapping like that for the STEMI samples as well
stemi_mapping <- create_anonymized_mapping(stemi@meta.data[stemi@meta.data[['orig.ident']] %in% c('stemi_v2', 'stemi_v3'), ], source_column = 'assignment.final', target_column = 'sample', target_prepend = 'S')
# merge the ll and stemi mapping
both_mapping <- rbind(ll_mappping, data.frame(sample = stemi_mapping[['sample']], geno_id = stemi_mapping[['assignment.final']]))
# replace or remove columns with the assignment
stemi@meta.data[['assignment.final']] <- both_mapping[match(stemi@meta.data[['assignment.final']], both_mapping[['geno_id']]), 'sample']
stemi@meta.data[['SNG.1ST']] <- both_mapping[match(stemi@meta.data[['SNG.1ST']], both_mapping[['geno_id']]), 'sample']
stemi@meta.data[['assignment_ll']] <- both_mapping[match(stemi@meta.data[['SNG.1ST']], both_mapping[['geno_id']]), 'sample']
stemi@meta.data[['BEST']] <- NULL

# metadata
metadata_gz <- gzfile(paste(counts_export_loc, 'metadata.tsv.gz', sep = ''))
write.table(cbind(data.frame(bc = rownames(stemi@meta.data)), stemi@meta.data), metadata_gz, row.names = F, col.names = T, quote = F, sep = '\t')

# subset our mapping to only assignments that are present
both_mapping <- both_mapping[both_mapping[['sample']] %in% stemi@meta.data[['assignment.final']], ]
# write the mapping
write.table(both_mapping, '/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/metadata/stemi_anon_id_to_geno_id.tsv', row.names = F, col.names = T)

# finally do just the raw data before filtering
stemi_v2_raw <- readRDS('/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/stemi_v2_raw_samples_20201110.rds')
stemi_v3_raw <- readRDS('/groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/stemi_v3_raw_samples_20201110.rds')
stemi_raw <- merge(stemi_v2_raw, stemi_v3_raw)
# features
features_gz_stemi_unfiltered <- gzfile(paste(counts_export_loc, 'stemi_unfiltered_features.tsv.gz', sep = ''))
write.table(data.frame(x = rownames(stemi_raw@assays$RNA@counts)), features_gz_stemi_unfiltered, row.names = F, col.names = F, quote = F)
# barcodes
barcodes_gz_stemi_unfiltered <- gzfile(paste(counts_export_loc, 'stemi_unfiltered_barcodes.tsv.gz', sep = ''))
write.table(data.frame(x = colnames(stemi_raw@assays$RNA@counts)), barcodes_gz_stemi_unfiltered, row.names = F, col.names = F, quote = F)
# and finally the count matrix
writeMM(stemi_raw@assays$RNA@counts, paste(counts_export_loc, 'stemi_unfiltered_matrix.mtx', sep = ''))
# matrix was manually zipped with bgzip /groups/umcg-franke-scrna/tmp02/releases/blokland-2020/v1/seurat/export/stemi_unfiltered_matrix.mtx
