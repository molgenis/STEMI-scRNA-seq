####################
# libraries        #
####################

library(MAST)
library(Seurat)
library(parallel)

####################
# Functions        #
####################

add_clinical_variables <- function(seurat_object, clinvar_object, clinical_vars_to_add, assignment_column_seurat='assignment.final', assignment_column_clinical='record_id', assignment_column_clinical_prepend='TEST_'){
  # subset the seurat object to just the metadata assignments
  meta.data <- seurat_object@meta.data[, c(assignment_column_seurat), drop = F]
  # make that a character type in case it was a factor
  meta.data[[assignment_column_seurat]] <- as.character(meta.data[[assignment_column_seurat]])
  # construct the assignment in the clinical variable object
  clinvar_object$assignment.seurat <- paste(assignment_column_clinical_prepend, as.character(clinvar_object[[assignment_column_clinical]]), sep = '')
  # subset the metadata to what we could possibly match, to save time
  meta.data <- meta.data[meta.data[[assignment_column_seurat]] %in% unique(clinvar_object$assignment.seurat), , drop = F]
  # grab the clinical variables in the order of the metadata by matching the ID column
  clinvar_matched <- clinvar_object[match(meta.data[[assignment_column_seurat]], clinvar_object$assignment.seurat), clinical_vars_to_add]
  # the order is now the same, so we can cbind
  meta.data <- cbind(meta.data, clinvar_matched)
  # now we can add it to the object, by subsetting the metadata by columns, it'll match by the rownames(cell barcodes)
  seurat_object <-AddMetaData(seurat_object, meta.data[, clinical_vars_to_add])
  return(seurat_object)
}

add_peak_ck_mb <- function(clinvar_object){
  # subset to the ckmb variables
  ckmb <- clinvar[, grep('ck_mb', colnames(clinvar_object))]
  # change the NA variables to -1
  ckmb[is.na(ckmb)] <- -1
  # get the maximum
  max_ckmb <- apply(ckmb, 1, max)
  # replace a -1 with NA, this means there was no real variable
  max_ckmb[max_ckmb == -1] <- NA
  # paste onto the existing object
  clinvar_object$peak_ck_mb <- max_ckmb
  return(clinvar_object)
}

seurat_object_to_singlecellassay <- function(seurat_object){
  # Convert Seurat object to SingleCellExperiment object; and then from SingleCellExperiment object to SingleCellAssay object.
  DefaultAssay(seurat_object) <- "RNA"
  sca <- as(as.SingleCellExperiment(seurat_object), 'SingleCellAssay')
  scaRaw <- sca
  # Filter out lowly variable genes
  sca <- sca[freq(sca)>0,] #returns the frequency of expression, i.e., the proportion of non-zero values in sc
  colData(sca)$wellKey<-colnames(sca)
  rowData(sca)$primerid <- rownames(rowData(sca))
  colData(sca) <- droplevels(colData(sca))
  colData(sca)$lane <- factor(colData(sca)$lane) # convert to factor
  colData(sca)$age_scaled <- scale(colData(sca)$age) # rescale continous variable (Age)
  cdr2 <- colSums(assay(sca)>0) #assay(sca) is working on assays(sca)$logcounts
  colData(sca)$cngeneson <- scale(cdr2)
  return(sca)
}

de_glmer.func <- function(sca_object, contrast, fixed_effects, random_effects, nagq0=F, simplified=T, alt_order=F){
  print(paste0('Testing: ', contrast))
  if(simplified){
    random_effects <- random_effects[random_effects!='lane']
  }
  random_effects.fmla <- paste(paste0('(1|',random_effects,')'),collapse='+') #try other nomenclature (=option 1, default)
  contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
  if(alt_order){
    contrast_fixed.fmla <- paste(c(fixed_effects,contrast),collapse='+')
  }
  zlm_vars <- paste0('~',paste(c(contrast_fixed.fmla,random_effects.fmla), collapse='+'))
  zlm_formula <- as.formula(zlm_vars)
  print(paste0('Fitting glmer: ',zlm_vars))
  if(nagq0){
    print(paste0('Trying to mitigate model failing to converge for some genes by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))'))
    zlmCond <- zlm(zlm_formula,
                   sca_object,
                   method='glmer', ebayes=FALSE, fitArgsD=list(nAGQ=0))
    summaryCond <- summary(zlmCond, doLRT=contrast, fitArgsD=list(nAGQ=0))
  }else{
    zlmCond <- zlm(zlm_formula,
                   sca_object,
                   method='glmer', ebayes=FALSE)
    summaryCond <- summary(zlmCond, doLRT=contrast)
  }
  summaryDt <- summaryCond$datatable
  z1_list <- list(zlm=zlmCond, dt=summaryDT)
  ## residuals -> not working with glmer!
  fixed_effects.fmla <- paste0('~',paste(fixed_effects,collapse='+'))
  print(paste0('Collecting residuals (with GLM): ',fixed_effects.fmla))
  fixed_effects.formula <- as.formula(fixed_effects.fmla)
  # window <- function(x1) lapply(assays(x1), function(x2) x2[, 1:2])
  ### total residuals of the response -> not working with glmer!
  z2 <- zlm(fixed_effects.formula, sca_object, hook=discrete_residuals_hook)
  z2_residuals <- collectResiduals(z2, sca_object)
  z2_summaryDt <- summary(z2, doLRT=contrast)$datatable
  z2_list <- list(zlm=z2, dt=z2_summaryDt, residuals=z2_residuals)
  # window(z1_residuals)
  z3 <- zlm(fixed_effects.formula, sca_object, hook=continuous_residuals_hook)
  z3_residuals <- collectResiduals(z3, sca_object)
  z3_summaryDt <- summary(z3, doLRT=contrast)$datatable
  z3_list <- list(zlm=z3, dt=z3_summaryDt, residuals=z3_residuals)
  # window(z2_residuals)
  z4 <- zlm(fixed_effects.formula, sca_object, hook=combined_residuals_hook)
  z4_residuals <- collectResiduals(z4, sca_object)
  z4_summaryDt <- summary(z4, doLRT=contrast)$datatable
  z4_list <- list(zlm=z4, dt=z4_summaryDt, residuals=z4_residuals)
  # window(z3_residuals)
  
  #residualsList <- list(discrete = assays(z1_residuals)$Residuals,
  #                      continuous = assays(z2_residuals)$Residuals,
  #                      combined = assays(z3_residuals)$Residuals)
  #summary_bygene <- lapply(residualsList, function(x) summary(t(x))) #check summary of each of the residuals (discrete, continuous and combined) per gene
  
  # summary_bygene
  ### partial residuals -> not working neither with glm
  contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
  contrast_fixed.formula <- as.formula(paste0('~',contrast_fixed.fmla))
  z5 <- zlm(contrast_fixed.formula, sca_object)
  partialScore(z5, 'MTL.value')
  z5_summaryDt <- summary(z5, doLRT=contrast)$datatable
  z5_list <- list(zlm=z5, dt=z5_summaryDt)
  ## save residuals
  #residuals.fn <- paste0(ct.dir, 'residuals_glm.rds')
  
  #print(paste0('Saving residuals with MAST GLM: ', residuals.fn))
  #saveRDS(residualsList, residuals.fn)
  
  # put it all in the list
  z_list <- list(z1=z1_list, z2=z2_list, z3=z3_list, z4=z4_list, z5=z5_list)
  
  return(z_list)
}

perform_MAST <- function(sca, contrast.vars=list('peak_ck_mb'), fixed_effects=c('cngeneson','gender','age', 'ck_mb'), random_effects=c('assignment.final','lane'), nagq0=F){
  # apply for each contrast
  de_glmer.list <- sapply(contrast.vars, function(i) de_glmer.func(sca_object = sca,
                                                                   contrast = i,
                                                                   fixed_effects = fixed_effects,
                                                                   random_effects = random_effects,
                                                                   nagq0 = nagq0), simplify = FALSE)
  # create result list
  fcHurdleSigs <- list
  # retrieve DEGs
  for(contrast in names(de_glmer.list)){
    summaryDt <- de_glmer.list[[contrast]]
    fcHurdle <- merge(summaryDt[contrast==contrast & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                      summaryDt[contrast==contrast & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], #logFC coefficients
                      by = 'primerid')
    fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    fcHurdleSig <- merge(fcHurdle[fdr<.05],
                         as.data.table(mcols(sca)),
                         by='primerid')
    # add to list
    fcHurdleSigs[[contrast]] <- fcHurdleSig
  }
  return(fcHurdleSigs)
}

glmer_per_condition <- function(seurat_object, conditions=c('Baseline', 't24h', 't8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), condition.column='timepoint.final', cell.type.column='cell_type_lowerres', contrast.vars=list('peak_ck_mb'), fixed_effects=c('cngeneson','gender','age', 'ck_mb'), random_effects=c('assignment.final','lane')){
  # setup result per condition
  result_per_condition <- list()
  # loop the conditions
  for(condition in intersect(conditions, unique(seurat_object@meta.data[[condition.column]]))){
    # subset to the condition
    seurat_object_condition <- seurat_object[, seurat_object@meta.data[[condition.column]] == condition]
    # onto the next function
    results <- glmer_per_cell_type(seurat_object_condition, cell_types=cell_types, cell.type.column=cell.type.column, contrast.vars=contrast.vars, fixed_effects=fixed_effects, random_effects=random_effects)
    # store the results
    result_per_condition[[condition]] <- results
  }
  return(result_per_condition)
}

glmer_per_condition_parallel <- function(seurat_object, conditions=c('Baseline', 't24h', 't8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), condition.column='timepoint.final', cell.type.column='cell_type_lowerres', contrast.vars=list('peak_ck_mb'), fixed_effects=c('cngeneson','gender','age', 'ck_mb'), random_effects=c('assignment.final','lane')){
  # set the conditions to use
  conditions_to_test <- intersect(conditions, unique(seurat_object@meta.data[[condition.column]]))
  # do a multithreaded loop over the list
  result_per_condition <- mclapply(conditions_to_test, function(x){
    # subset to the condition
    seurat_object_condition <- seurat_object[, seurat_object@meta.data[[condition.column]] == x]
    # onto the next function
    results <- glmer_per_cell_type_parallel(seurat_object_condition, cell_types=cell_types, cell.type.column=cell.type.column, contrast.vars=contrast.vars, fixed_effects=fixed_effects, random_effects=random_effects)
    return(results)
  })
  # set the conditions as the keys
  names(result_per_condition) <- conditions_to_test
  return(result_per_condition)
}


glmer_per_cell_type <- function(seurat_object, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell.type.column='cell_type_lowerres', contrast.vars=list('peak_ck_mb'), fixed_effects=c('cngeneson','gender','age', 'ck_mb'), random_effects=c('assignment.final','lane')){
  # setup result per cell_type
  result_per_cell_type <- list()
  # loop the cell types
  for(cell_type in intersect(cell_types, unique(seurat_object@meta.data[[cell.type.column]]))){
    # subset to the cell type
    seurat_object_cell_type <- seurat_object[, seurat_object@meta.data[[cell.type.column]] == cell_type]
    # convert to sca
    sca_ct <- seurat_object_to_singlecellassay(seurat_object_cell_type)
    
    # filter
    freq_expressed <- 0.1 #setting a minimum gene expression threshold based on genes that are found in at least 'freq_expressed' of the cells (proportion of non-zero cells)
    expressed_genes <- freq(sca_ct) > freq_expressed
    sca_ct <- sca_ct[expressed_genes,]
    
    # do the analysis
    mast_result <- perform_MAST(sca_ct, contrast.vars=contrast.vars, fixed_effects=fixed_effects, random_effects=random_effects)
    # store the result
    result_per_cell_type[[cell_type]] <- mast_result
  }
  return(result_per_cell_type)
}


glmer_per_cell_type_parallel <- function(seurat_object, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), cell.type.column='cell_type_lowerres', contrast.vars=list('peak_ck_mb'), fixed_effects=c('cngeneson','gender','age', 'ck_mb'), random_effects=c('assignment.final','lane')){
  # setup result per cell_type
  cell_types_to_test <- intersect(cell_types, unique(seurat_object@meta.data[[cell.type.column]]))
  # do a multithreaded loop over the list
  result_per_cell_type <- mclapply(cell_types_to_test, function(x){
    # subset to the cell type
    seurat_object_cell_type <- seurat_object[, seurat_object@meta.data[[cell.type.column]] == x]
    # convert to sca
    sca_ct <- seurat_object_to_singlecellassay(seurat_object_cell_type)
    
    # filter
    freq_expressed <- 0.1 #setting a minimum gene expression threshold based on genes that are found in at least 'freq_expressed' of the cells (proportion of non-zero cells)
    expressed_genes <- freq(sca_ct) > freq_expressed
    sca_ct <- sca_ct[expressed_genes,]
    
    # do the analysis
    mast_result <- perform_MAST(sca_ct, contrast.vars=contrast.vars, fixed_effects=fixed_effects, random_effects=random_effects)
    return(mast_result)
  })
  # set the cell_types as the keys
  names(result_per_cell_type) <- cell_types_to_test
  return(result_per_cell_type)
}

####################
# Main Code        #
####################

# storage spots
storage_read_part <- 'tmp01'
storage_write_part <- 'tmp01'

# object locations
object_loc <- paste('/groups/umcg-wijmenga/', storage_read_part, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/', sep = '')
cardio_object_loc <- paste(object_loc, 'cardio.integrated.20210301.rds', sep = '')

# location of the clinical variables
clinvar_loc <- paste('/groups/umcg-wijmenga/', storage_read_part, '/projects/1M_cells_scRNAseq/ongoing/Cardiology/clinical_parameters/scRNAseq_clinical_data_20200910_final_notextfields.tsv', sep = '')

# read the object
cardio.integrated <- readRDS(cardio_object_loc)

# read the clinical variables
clinvar <- read.table(clinvar_loc, sep = '\t', header=T, dec = ',')
# add the peak ckmb
clinvar <- add_peak_ck_mb(clinvar)
# add the peak ckmb and ckmb at t0
cardio.integrated <- add_clinical_variables(cardio.integrated, clinvar, clinical_vars_to_add = c('ck_mb', 'peak_ck_mb'))

# specifically for the monocyte
options(mc.cores=12)
de_ck_mb_monocyte <- glmer_per_condition_parallel(cardio.integrated, cell_types = c('monocyte'))
