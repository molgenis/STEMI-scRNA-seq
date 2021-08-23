#
# some sort of description of the file
#

#
# libraries used
#

library(Seurat)
library(Matrix)

#
# functions
#

# do doublet detection and removal
remove_doublets <- function(seurat_object, detection_method="souporcell"){
  if(detection_method == "souporcell"){
    # the status keyword states whether this was a singlet or not
    seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data['status']!= "singlet", "possible_doublet")
    # then select by this new row
    seurat_object <- subset(seurat_object, subset = possible_doublet == F)
  }
  else if(detection_method == "demuxlet"){
    # include based on singlet likelyhood
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1 < 25)
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
  }
  return(seurat_object)
}

do_pbmc_reference_mapping <- function(seurat_object, reference_object){
  # set to the assay used in the reference
  DefaultAssay(seurat_object) <- 'SCT'
  
  # find transfer anchors between the reference and the query, the query is your dataset
  anchors.seurat_object <- FindTransferAnchors(
    reference = reference_object,
    query = seurat_object,
    normalization.method = "SCT",
    reference.reduction = "spca",
    dims = 1:50,
    recompute.residuals = FALSE
  )
  
  # map the query to the reference, this will project your dataset onto the reference, and copy the labels with a certain amount of certainty (the score)
  seurat_object <- MapQuery(
    anchorset = anchors.seurat_object,
    query = seurat_object,
    reference = reference_object,
    refdata = list(
      celltype.l1 = "celltype.l1",
      celltype.l2 = "celltype.l2",
      predicted_ADT = "ADT"
    ),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
  return(seurat_object)
}


# creates a files per cell type, with in the columns the participants and in the rows the genes, with each cell the mean expression of that participant for that gene
create_features_files <- function(seurat_object, output_loc, cell_types_to_output = NULL, symbols.to.ensg=F, symbols.to.ensg.mapping = "genes.tsv", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", prepend_1 = F, metaqtl_format=F, make_plink_compat = F){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # get the individuals in the object
  individuals <- unique(seurat_object@meta.data[,sample.id.column.name])
  # remove any na individuals (sometimes caused if using demux assignments)
  individuals <- individuals[!is.na(individuals)]
  # if the user did not specify the cell types to output, just do all of them
  if(is.null(cell_types_to_output)){
    cell_types_to_output = unique(seurat_object@meta.data[[cell_type_column]])
  }
  # we will want to report on the number of cells per cell type, so that needs to be stored in a matrix
  cell_counts <- matrix(nrow=length(cell_types_to_output), ncol = length(individuals), 
                        dimnames = list(cell_types_to_output, individuals))
  # go through the cell types we want to output
  for (celltype in cell_types_to_output) {
    print(paste("doing", celltype))
    # grab the cells specific to that cell type
    cells_cell_type <- seurat_object[,seurat_object@meta.data[cell_type_column] == celltype]
    # create a matrix where we will store the mean expressions
    mean_expression_matrix <- matrix(nrow=nrow(cells_cell_type), ncol = length(individuals), 
                                     dimnames = list(rownames(cells_cell_type), individuals))
    # go through the individuals
    for (individual in individuals) {
      # take into account that there might be zero expression
      if (sum(is.na(cells_cell_type@meta.data[,sample.id.column.name]) | cells_cell_type@meta.data[,sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[,individual] <- 0
        cell_counts[celltype,individual] <- 0
      }
      else {
        # grab the cells of this cell type and this individual
        cells_cell_type_individual <- cells_cell_type[,cells_cell_type@meta.data[,sample.id.column.name] == individual]
        # calculate the mean expression for each gene
        mean_expression_matrix[,individual] <- rowMeans(cells_cell_type_individual)
        # get the number of cells of this cell type for this individual
        cell_counts[celltype,individual] <- ncol(cells_cell_type_individual)
      }
    }
    # remap symbols to ensg numbers if necessary
    if (symbols.to.ensg) {
      genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
      genes$V2 <- gsub("_", "-", make.unique(genes$V2))
      rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
    }
    if(prepend_1){
      # prepend the '1_' to the id 
      colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
    }
    if(make_plink_compat){
      # replace the underscores with dot if using plink for mapping
      colnames(mean_expression_matrix) <- gsub('_', '.', colnames(mean_expression_matrix))
    }
    
    if(metaqtl_format){
      gene_names <- data.frame(id=rownames(mean_expression_matrix))
      mean_expression_matrix <- cbind(gene_names, data.frame(mean_expression_matrix))
      # write our table of means for this cell type
      write.table(mean_expression_matrix, 
                  file = paste0(output_loc, celltype, "_expression", ".tsv"),
                  quote = F, sep = "\t", col.names = T, row.names = F)
    }
    else{
      # write our table of means for this cell type
      write.table(mean_expression_matrix, 
                  file = paste0(output_loc, celltype, "_expression", ".tsv"),
                  quote = F, sep = "\t", col.names = NA)
    }
    print(paste("finished", celltype))
  }
  # write out table of 
  write.table(cell_counts, 
              file = paste0(output_loc, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

# creates a files per cell type, with in the columns the participants and in the rows the genes, with each cell the mean expression of that participant for that gene, in each folder per condition
create_feature_files_per_condition <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", symbols.to.ensg = F, symbols.to.ensg.mapping="genes.tsv", prepend_1 = F, metaqtl_format=F, make_plink_compat=F){
  # subset to what we can actually use
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & !is.na(seurat_object@meta.data[[condition.column.name]]) & !is.na(seurat_object@meta.data[[sample.id.column.name]])]
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # unless we wnat only specific conditions, we'll just to those
  if(!(is.null(conditions_to_output))){
    conditions <- conditions_to_output
  }
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    print(paste("grabbing condition", condition))
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    print(paste("finished grabbing condition", condition))
    # do the feature file creation for this condition
    create_features_files(cells_condition, cell_types_to_output = cell_types_to_output, output_loc = output_loc_condition, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping=symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = prepend_1, metaqtl_format = metaqtl_format, make_plink_compat=make_plink_compat)
  }
}

create_feature_files_per_condition_combined <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output_a=c('UT'), conditions_to_output_b=c('Baseline', 't24h', 't8w'), condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", symbols.to.ensg = F, symbols.to.ensg.mapping="genes.tsv", prepend_1_a = T, prepend_1_b = F, metaqtl_format=F, make_plink_compat=F){
  # subset to what we can actually use
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & !is.na(seurat_object@meta.data[[condition.column.name]]) & !is.na(seurat_object@meta.data[[sample.id.column.name]])]
  # combined each A with each B
  for(condition_a in conditions_to_output_a){
    # with condition b
    for(condition_b in conditions_to_output_b){
      # subset to both
      seurat_object_conditions <- seurat_object[, seurat_object@meta.data[[condition.column.name]] %in% c(condition_a, condition_b)]
      # do magic on the participant column
      if(prepend_1_a){
        seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_a, sample.id.column.name] <- paste('1_', seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_a, sample.id.column.name], sep = '')
      }
      if(prepend_1_b){
        seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_b, sample.id.column.name] <- paste('1_', seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_b, sample.id.column.name], sep = '')
      }
      # paste the output location
      output_loc_full <- paste(output_loc, condition_a, '_', condition_b, '/', sep = '')
      # do the feature file creation for these conditions
      create_features_files(seurat_object_conditions, cell_types_to_output = cell_types_to_output, output_loc = output_loc_full, symbols.to.ensg=symbols.to.ensg, symbols.to.ensg.mapping=symbols.to.ensg.mapping, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = F, metaqtl_format = metaqtl_format, make_plink_compat = make_plink_compat)
    }
  }
}


# creates a files per cell type, with in the columns the participants and in the rows the genes, with each cell the mean expression of that participant for that gene, in each folder per condition
create_metadata_files_per_condition <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output = NULL, condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", prepend_1 = F){
  # grab the conditions
  conditions <- unique(seurat_object@meta.data[[condition.column.name]])
  # remove any na conditions (sometimes caused if using demux assignments)
  conditions <- conditions[!is.na(conditions)]
  # unless we wnat only specific conditions, we'll just to those
  if(!(is.null(conditions_to_output))){
    conditions <- conditions_to_output
  }
  # go through the conditions
  for(condition in conditions){
    # there should be a subdirectory per condition
    output_loc_condition <- paste0(output_loc, condition,"/")
    # there might be a prepended 'X' that we might want to remove
    if(startsWith(condition, "X")){
      output_loc_condition <- paste0(output_loc, substr(condition, 2, nchar(condition)),"/")
    }
    # grab only the cells with this condition
    print(paste("grabbing condition", condition))
    cells_condition <- seurat_object[,seurat_object@meta.data[condition.column.name] == condition]
    print(paste("finished grabbing condition", condition))
    # do the feature file creation for this condition
    create_metadata_files(cells_condition, cell_types_to_output = cell_types_to_output, output_loc = output_loc_condition, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = prepend_1)
  }
}


create_metadata_files_per_condition_combined <- function(seurat_object, output_loc, cell_types_to_output = NULL, conditions_to_output_a=c('UT'), conditions_to_output_b=c('Baseline', 't24h', 't8w'), condition.column.name = "timepoint.final", sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", symbols.to.ensg = F, symbols.to.ensg.mapping="genes.tsv", prepend_1_a = T, prepend_1_b = F, metaqtl_format=F, make_plink_compat=F){
  # subset to what we can actually use
  seurat_object <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & !is.na(seurat_object@meta.data[[condition.column.name]]) & !is.na(seurat_object@meta.data[[sample.id.column.name]])]
  # combined each A with each B
  for(condition_a in conditions_to_output_a){
    # with condition b
    for(condition_b in conditions_to_output_b){
      # subset to both
      seurat_object_conditions <- seurat_object[, seurat_object@meta.data[[condition.column.name]] %in% c(condition_a, condition_b)]
      # do magic on the participant column
      if(prepend_1_a){
        seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_a, sample.id.column.name] <- paste('1_', seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_a, sample.id.column.name], sep = '')
      }
      if(prepend_1_b){
        seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_b, sample.id.column.name] <- paste('1_', seurat_object_conditions@meta.data[seurat_object_conditions@meta.data[[condition.column.name]] == condition_b, sample.id.column.name], sep = '')
      }
      # paste the output location
      output_loc_full <- paste(output_loc, condition_a, '_', condition_b, '/', sep = '')
      # do the feature file creation for this condition
      create_metadata_files(seurat_object_conditions, cell_types_to_output = cell_types_to_output, output_loc = output_loc_full, sample.id.column.name=sample.id.column.name, cell_type_column=cell_type_column, assay=assay, prepend_1 = F)
    }
  }
}


create_metadata_files <- function(seurat_object,  output_loc, cell_types_to_output = NULL, sample.id.column.name="assignment.final", cell_type_column = "cell_type", assay = "RNA", prepend_1 = F){
  # set the assay
  DefaultAssay(seurat_object) <- assay
  # if the user did not specify the cell types to output, just do all of them
  if(is.null(cell_types_to_output)){
    cell_types_to_output = unique(seurat_object@meta.data[[cell_type_column]])
  }
  # go through the cell types we want to output
  for (celltype in cell_types_to_output) {
    # get the individuals in the object
    individuals <- unique(seurat_object@meta.data[seurat_object@meta.data[[cell_type_column]] == celltype, sample.id.column.name])
    # remove any na individuals (sometimes caused if using demux assignments)
    individuals <- individuals[!is.na(individuals)]
    # create a matrix where we will store the data
    metadata_matrix <- matrix(nrow=7, ncol = length(individuals), 
                              dimnames = list(c('chem', 'lane', 'ncells', 'avg_numi', 'gender', 'age', 'status'), individuals))
    # go through the individuals
    for (individual in individuals) {
      # get the data for the individual
      try({
        # subset to cells of individual
        seurat_object_participant <- seurat_object[, seurat_object@meta.data[[sample.id.column.name]] == individual &
                                                     seurat_object@meta.data[[cell_type_column]] == celltype]
        # get static data
        chem <- unique(seurat_object_participant@meta.data$chem)[1]
        lane <- unique(seurat_object_participant@meta.data$lane)[1]
        gender <- unique(seurat_object_participant@meta.data$gender)[1]
        age <- unique(seurat_object_participant@meta.data$age)[1]
        ncells <- nrow(seurat_object_participant@meta.data)
        status <- unique(seurat_object_participant@meta.data[['timepoint.final']])[1]
        # get the average number of umis
        avg_umi <- NULL
        if(assay == 'SCT'){
          avg_umi <- mean(seurat_object_participant@meta.data$nCount_SCT)
        }
        else if(assay == 'RNA'){
          avg_umi <- mean(seurat_object_participant@meta.data$nCount_RNA)
        }
        else{
          print('unknown assay, using RNA')
          avg_umi <- mean(seurat_object_participant@meta.data$nCount_RNA)
        }
        # setting data
        metadata_matrix['chem', individual] <- chem
        metadata_matrix['lane', individual] <- lane
        metadata_matrix['gender', individual] <- gender
        metadata_matrix['age', individual] <- age
        metadata_matrix['ncells', individual] <- ncells
        metadata_matrix['avg_numi', individual] <- avg_umi
        metadata_matrix['status', individual] <- status
      })
    }
    # write the results
    if(prepend_1){
      # prepend the '1_' to the id 
      colnames(metadata_matrix) <- paste0("1_", colnames(metadata_matrix))
    }
    # add the rownames as a variable
    var_names <- data.frame(id=c('chem', 'lane', 'ncells', 'avg_numi', 'gender', 'age', 'status'))
    metadata_matrix <- cbind(var_names, data.frame(metadata_matrix))
    # write our table of means for this cell type
    write.table(metadata_matrix, 
                file = paste0(output_loc, gsub(' ', '_', celltype), "_metadata", ".tsv"),
                quote = F, sep = "\t", col.names = T, row.names = F)
  }
}

metadata_file_to_coded_file <- function(metadata_file_loc, metadata_output_loc=NULL,vars_to_keep=NULL, vars_to_convert=NULL){
  # read the original file
  metadata <- read.table(metadata_file_loc, sep = '\t', header=T, stringsAsFactors = F)
  # create new metadata
  metadata_new <- NULL
  # check if the user gave us any commands
  variables <- vars_to_keep
  if(is.null(variables)){
    variables <- metadata$id
  }
  varsconv <- vars_to_convert
  if(is.null(varsconv)){
    varsconv <- metadata$id
  }
  # here's what we'll keep
  metadata_new <- metadata[metadata$id %in% variables, , drop=F]
  # and now do some replacing
  for(variable_to_convert in vars_to_convert){
    # get the variables in the row
    variable_row <- as.character(as.vector(unlist(metadata[metadata$id == variable_to_convert, 2:ncol(metadata)])))
    # check the unique variables
    variable_unique <- unique(variable_row)
    # turn binary, skipping the first one, because that would cause duplicate describing rows
    if(length(variable_unique) > 1){
      for(i in 2:length(variable_unique)){
        # check the option
        variable <- as.character(variable_unique[i])
        # initiate a row, where nothing is this specific variable
        binary_row <- rep(0, times = length(variable_row))
        # then set to 1 where in the original data, this position held this variable
        binary_row[variable_row == variable] <- 1
        # add as a binary row
        #metadata_new[paste(variable_to_convert, variable, sep = '_'), ] <- c(binary_row)
        row <- c(paste(variable_to_convert, variable, sep = '_'), binary_row)
        metadata_new <- rbind(metadata_new, row)
      }
      # remove the original row 
      metadata_new <- metadata_new[metadata_new$id != variable_to_convert, ]
    }
    else{
      # if the variable is the same, add it with all values being '1', so true
      row <- c(paste(variable_to_convert, variable_unique, sep = '_'), rep(1, times = length(variable_row)))
      metadata_new <- rbind(metadata_new, row)
      # remove the original row 
      metadata_new <- metadata_new[metadata_new$id != variable_to_convert, ]
    }
  }
  # write the result
  output_location <- metadata_output_loc
  if(is.null(output_location)){
    output_location <- paste(gsub('(tsv$)|(txt$)', '', metadata_file_loc), paste(vars_to_keep, collapse='_', sep = ''), '.tsv', sep='')
  }
  print(metadata_new)
  #metadata$id <- vars_to_keep
  # write our table of means for this cell type
  write.table(metadata_new, 
              file = output_location,
              quote = F, sep = "\t", col.names = T, row.names = F)
}




#
# main code
#

# where to store the objects
object_loc <- "/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/"
cells_1M_loc <- paste(object_loc, '1M_cells_seurat_object_hg19_soup_raw.rds', sep = '')
cardio.integrated_loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')
v2_combined_loc <- paste(object_loc, 'combined.v2.20210629.rds', sep = '')
v3_combined_loc <- paste(object_loc, 'combined.v3.20210629.rds', sep = '')

# for the inhouse eQTL-mapping pipeline, we currently have ENSG numbers instead of gene numbers
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"


# we have done the strategy of adding idents etc in stemi_preprocess_HC.R, so we can just load the object
cells_1M <- readRDS(cells_1M_loc)
# remove these lanes, they were stimulated (original script is fixed, but just to be safe)
cells_1M <- cells_1M[, !(cells_1M@meta.data$lane %in% c("190123_lane1", "190123_lane2"))]
# remove the doublets
cells_1M <- remove_doublets(cells_1M)
# grab just UT
cells_1M <- subset(cells_1M, subset = timepoint.final == 'UT')
# split into two objects
cells_1M_v2 <- subset(cells_1M, subset = chem == 'V2')
cells_1M_v3 <- subset(cells_1M, subset = chem == 'V3')
# clear some memory
rm(cells_1M)
cells_1M_v2[["percent.mt"]] <- PercentageFeatureSet(cells_1M_v2, pattern = "^MT-")
cells_1M_v3[["percent.mt"]] <- PercentageFeatureSet(cells_1M_v3, pattern = "^MT-")
# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
cells_1M_v2 <- subset(cells_1M_v2, subset = nFeature_RNA > 200 & percent.mt < 8 & HBB < 10)
cells_1M_v3 <- subset(cells_1M_v3, subset = nFeature_RNA > 200 & percent.mt < 15 & HBB < 10)

# read the integrated object
cardio.integrated <- readRDS(cardio.integrated_loc)
# subset to the STEMI datasets
stemi_v2 <- subset(cardio.integrated, subset = orig.ident == 'stemi_v2')
stemi_v3 <- subset(cardio.integrated, subset = orig.ident == 'stemi_v3')
# clear up some memory
rm(cardio.integrated)



# add Seurat objects together
v2_combined <- merge(cells_1M_v2, stemi_v2)
v3_combined <- merge(cells_1M_v3, stemi_v3)
# clear up some memory
rm(cells_1M_v2)
rm(cells_1M_v3)
rm(stemi_v2)
rm(stemi_v3)

# now we have the complete object per chemistry, we should normalize them together
#v2_combined <- SCTransform(v2_combined, vars.to.regress = c('percent.mt'))
#v3_combined <- SCTransform(v3_combined, vars.to.regress = c('percent.mt'))
v2_combined <- SCTransform(v2_combined)
# do old style normalization as well
DefaultAssay(v2_combined) <- 'RNA'
v2_combined <- NormalizeData(v2_combined)
# back to SCT, which we will actually use
DefaultAssay(v2_combined) <- 'SCT'
# standard workflow stuff
v2_combined <- RunPCA(v2_combined)
v2_combined <- RunUMAP(v2_combined, dims=1:30, reduction='pca')
v2_combined <- FindNeighbors(v2_combined)
v2_combined <- FindClusters(v2_combined, resolution=1.2)
saveRDS(v2_combined, v2_combined_loc)
# for both sets
v3_combined <- SCTransform(v3_combined)
v3_combined <- RunPCA(v3_combined)
v3_combined <- RunUMAP(v3_combined, dims=1:30, reduction='pca')
v3_combined <- FindNeighbors(v3_combined)
v3_combined <- FindClusters(v3_combined, resolution=1.2)
saveRDS(v3_combined, v3_combined_loc)

# we have Azimuth predictions for the 1M UT cells, that we calculated before
v3_predictions_ut <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/v3_ut_20201106_azimuth_predictions.tsv', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
v2_predictions_ut <- read.table('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell-type-classifying/v2_ut_20201029_azimuth_predictions.tsv', sep = '\t', header = T, row.names = 1, stringsAsFactors = F)
for(column in colnames(v3_predictions_ut)){
  # add as a new column
  v3_combined <- AddMetaData(object = v3_combined, metadata = v3_predictions_ut[, column, drop = F], col.name = paste('new', column, sep='.'))
  # overwrite the existing column with the new column, where the new one is not empty (so only overwriting the 1M entries)
  v3_combined@meta.data[!is.na(v3_combined@meta.data[[paste('new', column, sep='.')]]), column] <- v3_combined@meta.data[!is.na(v3_combined@meta.data[[paste('new', column, sep='.')]]), paste('new', column, sep='.')]
  # remove that new column
  v3_combined@meta.data[[paste('new', column, sep='.')]] <- NULL
}
# for v2 as well
for(column in colnames(v2_predictions_ut)){
  v2_combined <- AddMetaData(object = v2_combined, metadata = v2_predictions_ut[, column, drop = F], col.name = paste('new', column, sep='.'))
  v2_combined@meta.data[!is.na(v2_combined@meta.data[[paste('new', column, sep='.')]]), column] <- v2_combined@meta.data[!is.na(v2_combined@meta.data[[paste('new', column, sep='.')]]), paste('new', column, sep='.')]
  v2_combined@meta.data[[paste('new', column, sep='.')]] <- NULL
}
# we'll convert it to a String, just to be sure that the next few steps are easier
v2_combined@meta.data$cell_type <- as.character(v2_combined@meta.data$cell_type)
v3_combined@meta.data$cell_type <- as.character(v3_combined@meta.data$cell_type)
# we will use the l2 predicted cell types (for the most part), for the classification of the 1M UT
v2_combined@meta.data[is.na(v2_combined@meta.data$cell_type), ]$cell_type <- as.character(v2_combined@meta.data[is.na(v2_combined@meta.data$cell_type), ]$predicted.celltype.l2)
v3_combined@meta.data[is.na(v3_combined@meta.data$cell_type), ]$cell_type <- as.character(v3_combined@meta.data[is.na(v3_combined@meta.data$cell_type), ]$predicted.celltype.l2)
# add another string to denote these as unclassified
v2_combined@meta.data[is.na(v2_combined@meta.data$cell_type), 'cell_type'] <- 'unclassified'
v3_combined@meta.data[is.na(v3_combined@meta.data$cell_type), 'cell_type'] <- 'unclassified'

# harmonize the cell types again
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
v2_combined@meta.data[v2_combined@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
# spaces in variables is inconvenient in a lot of places, so we'll replace these with underscores
v2_combined@meta.data$cell_type <- gsub(' ', '_', v2_combined@meta.data$cell_type)
# for v3 as well
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'NK', ]$cell_type <- 'NKdim'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'NK_CD56bright', ]$cell_type <- 'NKbright'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'CD14 Mono', ]$cell_type <- 'cMono'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'CD16 Mono', ]$cell_type <- 'ncMono'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'Plasmablast', ]$cell_type <- 'plasmablast'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'Platelet', ]$cell_type <- 'platelet'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'Eryth', ]$cell_type <- 'eryth'
v3_combined@meta.data[v3_combined@meta.data$cell_type == 'Doublet', ]$cell_type <- 'doublet'
# spaces in variables is inconvenient in a lot of places, so we'll replace these with underscores
v3_combined@meta.data$cell_type <- gsub(' ', '_', v3_combined@meta.data$cell_type)

# we will define a lower resolution cell type as well, we need to create some groupings for this
cd4t <- c('Treg', 'CD4_Naive', 'CD4_TCM', 'CD4_TEM', 'CD4_CTL', 'CD4_Proliferating')
cd8t <- c('MAIT', 'CD8_Naive', 'CD8_TCM', 'CD8_TEM', 'CD8_Proliferating')
t_other <- c('dnT', 'gdT', 'ILC')
nk <- c('NKdim', 'NKbright', 'NK_Proliferating')
monocyte <- c('cMono', 'ncMono')
dc <- c('cDC1', 'cDC2', 'pDC', 'ASDC')
b <- c('B_naive', 'B_intermediate', 'B_memory')
# add the new column by copying the higher res first, 
v2_combined@meta.data$cell_type_lowerres <- v2_combined@meta.data$cell_type
# in this new column, overwrite them to have the lower resolution
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% t_other, ]$cell_type_lowerres <- 'T_other'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% nk, ]$cell_type_lowerres <- 'NK'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% dc, ]$cell_type_lowerres <- 'DC'
v2_combined@meta.data[v2_combined@meta.data$cell_type %in% b, ]$cell_type_lowerres <- 'B'

# add the new column by copying the higher res first, 
v3_combined@meta.data$cell_type_lowerres <- v3_combined@meta.data$cell_type
# in this new column, overwrite them to have the lower resolution
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% cd4t, ]$cell_type_lowerres <- 'CD4T'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% cd8t, ]$cell_type_lowerres <- 'CD8T'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% t_other, ]$cell_type_lowerres <- 'T_other'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% nk, ]$cell_type_lowerres <- 'NK'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% monocyte, ]$cell_type_lowerres <- 'monocyte'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% dc, ]$cell_type_lowerres <- 'DC'
v3_combined@meta.data[v3_combined@meta.data$cell_type %in% b, ]$cell_type_lowerres <- 'B'

saveRDS(v2_combined, paste(object_loc, 'combined.v2.20210629.ct.rds', sep = ''))
saveRDS(v3_combined, paste(object_loc, 'combined.v3.20210629.ct.rds', sep = ''))

# create the feature files that were combined
#create_feature_files_per_condition_combined(v2_combined, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v2/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
#create_feature_files_per_condition_combined(v3_combined, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v3/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition_combined(v2_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v2/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition_combined(v3_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v3/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
# and only the 1M UT ones
#create_feature_files_per_condition(v2_combined[, v2_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v2/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
#create_feature_files_per_condition(v3_combined[, v3_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v3/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition(v2_combined[, v2_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v2/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition(v3_combined[, v3_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_v3/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')

# now combine the objects and do the v2 and v3 at the same time, we will merge with data, as we can use the chem as covariate later
all_combined <- merge(v2_combined, v3_combined, merge.data = T)
#create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
#create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = T, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT')
# add a mock bulk column
all_combined@meta.data$bulk <- 'bulk'
create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_limix/', cell_types_to_output=c('bulk'), cell_type_column = 'bulk', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=T)
create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_limix/', cell_types_to_output=c('bulk'), cell_type_column = 'bulk', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=T)
create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_limix/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=T)
create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_limix/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=T)

# create metadata for MatrixEQTL
create_metadata_files_per_condition(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/', cell_type_column = 'cell_type_lowerres')
create_metadata_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/', cell_type_column = 'cell_type_lowerres')
# add a mock bulk column
all_combined@meta.data$bulk <- 'bulk'
create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/', cell_types_to_output=c('bulk'), cell_type_column = 'bulk', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=F, metaqtl_format=T)
create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/', cell_types_to_output=c('bulk'), cell_type_column = 'bulk', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=F, metaqtl_format=T)
create_feature_files_per_condition_combined(all_combined, '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1_a = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=F, metaqtl_format=T)
create_feature_files_per_condition(all_combined[, all_combined@meta.data$timepoint.final == 'UT'], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/', cell_types_to_output=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell_type_column = 'cell_type_lowerres', prepend_1 = F, symbols.to.ensg = T, symbols.to.ensg.mapping = gene_to_ens_mapping, assay = 'SCT', make_plink_compat=F, metaqtl_format=T)


for(condition in c('Baseline', 't24h', 't8w', 'UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w')){
  print(condition)
  for(cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
    print(cell_type)
    loc <- paste('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/', condition, '/', cell_type, '_metadata.tsv', sep = '')
    metadata_file_to_coded_file(loc, vars_to_keep = c('chem', 'status'), vars_to_convert = c('chem', 'status'))
  }
}
