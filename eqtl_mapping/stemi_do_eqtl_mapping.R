# load libraries
library(MatrixEQTL)
library(data.table)
#library(foreach)
#library(doMC)

# functions
do_QTL_mapping <- function(
  SNP_file_name, # Genotype file name
  expression_file_name, # Gene expression file name
  snpspos, # dataframe containing the snp positions
  genepos, # dataframe containing the gene positions
  output_file_name_cis, # Output file name
  covariates_file_name=character(), # Covariates file name
  output_file_name_tra=tempfile(), # Output file name
  pvOutputThreshold_cis=1, # Only associations significant at this level will be saved
  pvOutputThreshold_tra=0, # Only associations significant at this level will be saved
  useModel=modelLINEAR, # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
  cisDist=1e6, # Distance for local gene-SNP pairs
  errorCovariance=numeric(), # Error covariance matrix
  snps_fileDelimiter='\t',# the TAB character
  snps_fileOmitCharacters='NA', # denote missing values;
  snps_fileSkipRows=1, # one row of column labels
  snps_fileSkipColumns=1, # one column of row labels
  snps_fileSliceSize=2000, # read file in slices of 2,000 rows
  gene_fileDelimiter='\t',# the TAB character
  gene_fileOmitCharacters='NA', # denote missing values;
  gene_fileSkipRows=1, # one row of column labels
  gene_fileSkipColumns=1, # one column of row labels
  gene_fileSliceSize=2000, # read file in slices of 2,000 rows
  cvrt_fileDelimiter='\t',# the TAB character
  cvrt_fileOmitCharacters='NA', # denote missing values;
  cvrt_fileSkipRows=1, # one row of column labels
  cvrt_fileSkipColumns=1, # one column of row labels
  verbose=TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
){
  ## Load genotype data
  snps = SlicedData$new();
  snps$fileDelimiter = snps_fileDelimiter      # the TAB character
  snps$fileOmitCharacters = snps_fileOmitCharacters; # denote missing values;
  snps$fileSkipRows = snps_fileSkipRows;          # one row of column labels
  snps$fileSkipColumns = snps_fileSkipColumns;       # one column of row labels
  snps$fileSliceSize = snps_fileSliceSize;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  gene = SlicedData$new();
  gene$fileDelimiter = gene_fileDelimiter;      # the TAB character
  gene$fileOmitCharacters = gene_fileOmitCharacters; # denote missing values;
  gene$fileSkipRows = gene_fileSkipRows;          # one row of column labels
  gene$fileSkipColumns = gene_fileSkipColumns;       # one column of row labels
  gene$fileSliceSize = gene_fileSliceSize;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  ## Load covariates
  cvrt = SlicedData$new();
  cvrt$fileDelimiter = cvrt_fileDelimiter;      # the TAB character
  cvrt$fileOmitCharacters = cvrt_fileOmitCharacters; # denote missing values;
  cvrt$fileSkipRows = cvrt_fileSkipRows;          # one row of column labels
  cvrt$fileSkipColumns = cvrt_fileSkipColumns;       # one column of row labels
  if(length(covariates_file_name)>0) {
    cvrt$LoadFile(covariates_file_name);
  }
  # run the actual application
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name     = output_file_name_tra,
    pvOutputThreshold     = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = verbose,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = min.pv.by.genesnp,
    noFDRsaveMemory = noFDRsaveMemory);
  # release resources for the temporary file
}

determine_fdr <- function(output_file_name_cis, permutation_rounds, do_smallest_per_gene=T){
  # we will add the permutation rounds together
  permutation_table <- NULL
  # we'll store the SNP and gene info separately, that's easier for the rowmeans
  snp_gene <- NULL
  for(i in 1:permutation_rounds){
    # read the current round
    #permutation_round <- fread(paste(output_file_name_cis, '.permuted.', i, sep = ''), sep = '\t', header = T)
    # reset the order
    permutation_round <- read.table(paste(output_file_name_cis, '.permuted.', i, sep = ''), sep = '\t', header = T, stringsAsFactors = F)
    rownames(permutation_round) <- paste(permutation_round[['SNP']], permutation_round[['gene']], sep = '_')
    # paste onto the rest of the output
    if(is.null(permutation_table)){
      #permutation_table <- data.frame(perm1 = permutation_round[['p-value']])
      permutation_table <- permutation_round[, c('p.value'), drop = F]
      permutation_table[['perm1']] <- permutation_round[['p.value']]
      permutation_round[['p.value']] <- NULL
      # set the rownames here as well
      #rownames(permutation_table) <- rownames(permutation_round)
      # we'll need this one later
      snp_gene <- permutation_round[, c('SNP', 'gene')]
    }
    else{
      # grab the same snp-gene pairs
      permutation_table[[paste('perm', i, sep = '')]] <- permutation_round[rownames(permutation_table), c('p.value'), drop = F]
    }
  }
  # calculate the average p-value
  mean_p_permuted=rowMeans(permutation_table)
  # add back the snp and gene info
  permuted_table <- snp_gene
  permuted_table[['perm']] <- mean_p_permuted
  # order by significance
  permuted_table <- permuted_table[order(permuted_table[['perm']]), ]
  # filter to the best value per gene, if requested
  if(do_smallest_per_gene){
    # we ordered the table, so the first entry found for a gene, should always be the most significant one
    permuted_table <- permuted_table[match(unique(permuted_table[['gene']]), permuted_table[['gene']]), ]
  }
  # now read the non-permuted file
  #output_file <- fread(paste(output_file_name_cis), sep = '\t', header = T)
  output_file <- read.table(paste(output_file_name_cis), sep = '\t', header = T, stringsAsFactors = F)
  # go through each row
  output_file[['permuted_fdr']] <- apply(output_file, 1, function(x){
    # grab the p-value
    p_real <- x[['p.value']]
    # get the number of permuted p values
    nr_permuted_p <- nrow(permuted_table)
    # check how many of these p values are smaller than the actual p
    nr_permuted_p_smaller <- nrow(permuted_table[permuted_table[['perm']] < p_real, ])
    # calculate the fraction of permuted Ps, that are smaller than your actual p
    perm_fdr <- nr_permuted_p_smaller / nr_permuted_p
    # that fraction is the chance that your p could come from the permuted distribution
    return(perm_fdr)
  })
  # get the file without the extention
  output_file_name_cis_fdred <- sub('\\..[^\\.]*$', '', output_file_name_cis)
  # add new extention
  output_file_name_cis_fdred <- paste(output_file_name_cis_fdred, '.fdr.tsv', sep = '')
  # write the file
  write.table(output_file, output_file_name_cis_fdred, sep = '\t', row.names = F, col.names = T, quote = F)
}


run_qtl_mapping <- function(features_loc_ct_cond, output_file_name_cis_ct_cond, covariates_file_loc_ct_cond, snps, snpspos, genepos, maf=0.1, permutation_rounds = 0, permute_in_covar_group=NULL, do_smallest_per_gene=T, gene_confinement=NULL, snp_confinement=NULL){
  # read covariate data
  covariates_ct_cond <- fread(covariates_file_loc_ct_cond, sep = '\t', header = T, stringsAsFactors=FALSE)
  # read the expression data
  expressions_ct_cond <- fread(features_loc_ct_cond, sep = '\t', header=T, stringsAsFactors=FALSE)
  # subset to confined ones, if requested
  if(!is.null(gene_confinement)){
    expressions_ct_cond <- expressions_ct_cond[!is.na(expressions_ct_cond$id) & expressions_ct_cond$id %in% gene_confinement, ]
  }
  if(!is.null(snp_confinement)){
    snps <- snps[!is.na(snps$id) & snps$id %in% snp_confinement, ]
  }

  # get the participants that we have both expression and snps data for
  participants_ct_cond <- intersect(colnames(snps)[2:ncol(snps)], colnames(expressions_ct_cond)[2:ncol(expressions_ct_cond)])
  # get also overlap with the covariates data
  participants_ct_cond <- intersect(participants_ct_cond, colnames(covariates_ct_cond)[2:ncol(covariates_ct_cond)])
  # perform subsetting
  snps_ct_cond <- snps[, c('id', participants_ct_cond), with = F]
  expressions_ct_cond <- expressions_ct_cond[, c('id', participants_ct_cond), with = F]
  covariates_ct_cond <- covariates_ct_cond[, c('id', participants_ct_cond), with = F]

  # remove the covariates which show no variation
  covariates_ct_cond <- covariates_ct_cond[apply(covariates_ct_cond, 1, function(x){length(unique(x)) > 1}), ]
  
  # remove SNPs with missing values
  snps_ct_cond <- snps_ct_cond[rowSums(is.na(snps_ct_cond)) == 0 & rowSums(snps_ct_cond[, 2:ncol(snps_ct_cond)] < 0) == 0, ]
  
  # doublecheck SNPs
  all_env <<- snps_ct_cond
  
  # filter snps by maf, which of course can work both ways
  maf_reverse <- 1-maf
  snps_ct_cond <- snps_ct_cond[rowSums(snps_ct_cond[, 2:ncol(snps_ct_cond)])/(ncol(snps_ct_cond)-1)/2 >= maf & rowSums(snps_ct_cond[, 2:ncol(snps_ct_cond)])/(ncol(snps_ct_cond)-1)/2 <= maf_reverse, ]
  
  # now write the filtered files to temporary storage
  SNP_file_name_ct_cond <- tempfile()
  expression_file_name_ct_cond <- tempfile()
  covariates_file_name_ct_cond <- tempfile()
  print('writing the filtered files')
  write.table(snps_ct_cond, SNP_file_name_ct_cond, sep = '\t', row.names = F, col.names = T)
  write.table(expressions_ct_cond, expression_file_name_ct_cond, sep = '\t', row.names = F, col.names = T)
  write.table(covariates_ct_cond, covariates_file_name_ct_cond, sep = '\t', row.names = F, col.names = T)

  # do actual mapping
  print('start actual mapping')
  do_QTL_mapping(
    SNP_file_name=SNP_file_name_ct_cond, # Genotype file name
    expression_file_name=expression_file_name_ct_cond, # Gene expression file name
    snpspos=snpspos, # dataframe containing the snp positions
    genepos=genepos, # dataframe containing the gene positions
    output_file_name_cis=output_file_name_cis_ct_cond, # Output file name
    covariates_file_name=covariates_file_name_ct_cond, # Covariates file name
    cisDist = 100000
  )
  # do permuted mappings as well
  for(i in 1:permutation_rounds){
  #cl <- parallel::makeCluster(4)
  #doParallel::registerDoParallel(cl)
  #foreach(i=1:permutation_rounds) %dopar% {
    print(paste('starting permutation', i))
    # we will have a permuted expression file for each round
    expression_file_name_ct_cond_permuted <- tempfile()
    # and matched permuted covariates
    covariate_file_name_ct_cond_permuted <- tempfile()
    # and we will output somewhere
    output_file_name_cis_ct_cond_permuted <- paste(output_file_name_cis_ct_cond, '.permuted.', i, sep = '')
    # do label swapping per covariate group, if requested
    if(!is.null(permute_in_covar_group)){
      # start with the unpermuted file
      expressions_ct_cond_permuted <- expressions_ct_cond
      # covariates as well
      covariates_ct_cond_permuted <- covariates_ct_cond
      # go through each permutable category
      for(group in unique(covariates_ct_cond[[permute_in_covar_group]])){
        # get the samples with this group
        partipants_group <- covariates_ct_cond[covariates_ct_cond[[permute_in_covar_group]] == group, 1]
        # get the location of these participants in the expression file
        participant_locations <- match(partipants_group, colnames(covariates_ct_cond))
        # randomly shuffle these positions
        participant_locations_shuffled <- sample(x = participant_locations, size = length(participant_locations))
        # now replace the entries we have for this group, with the ones we sampled
        expressions_ct_cond_permuted[, participant_locations] <- expressions_ct_cond_permuted[, participant_locations_shuffled]
        # permute the covariates in the same way so they are still matched to the expression (actually permuting genotype)
        covariates_ct_cond_permuted[, participant_locations] <- covariates_ct_cond_permuted[, participant_locations_shuffled]
      }
      # and in the end, set the colnames like nothing happened
      colnames(expressions_ct_cond_permuted) <- colnames(expressions_ct_cond)
      colnames(covariates_ct_cond_permuted) <- colnames(covariates_ct_cond)
      # then write the file
      write.table(expressions_ct_cond_permuted, expression_file_name_ct_cond_permuted, sep = '\t', row.names = F, col.names = T)
      write.table(covariates_ct_cond_permuted, covariate_file_name_ct_cond_permuted, sep = '\t', row.names = F, col.names = T)
    }
    # do it all at once
    else{
      # grab randomly from the column names (all after the first), which should be the participants
      random_sampling <- c(1 , sample(x = 2:ncol(expressions_ct_cond), size = ncol(expressions_ct_cond) - 1))
      # grab participants 
      expressions_ct_cond_permuted <- expressions_ct_cond[, random_sampling]
      covariates_ct_cond_permuted <- covariates_ct_cond[, random_sampling]
      # set as if we did not change anything
      colnames(expressions_ct_cond_permuted) <- colnames(expressions_ct_cond)
      colnames(covariates_ct_cond_permuted) <- colnames(covariates_ct_cond)
      # write as file
      write.table(expressions_ct_cond_permuted, expression_file_name_ct_cond_permuted, sep = '\t', row.names = F, col.names = T)
      write.table(covariates_ct_cond_permuted, covariate_file_name_ct_cond_permuted, sep = '\t', row.names = F, col.names = T)
    }
    do_QTL_mapping(
      SNP_file_name=SNP_file_name_ct_cond, # Genotype file name
      expression_file_name=expression_file_name_ct_cond_permuted, # Gene expression file name
      snpspos=snpspos, # dataframe containing the snp positions
      genepos=genepos, # dataframe containing the gene positions
      output_file_name_cis=output_file_name_cis_ct_cond_permuted, # Output file name
      covariates_file_name=covariate_file_name_ct_cond_permuted, # Covariates file name
      cisDist = 100000
    )
  }
  #parallel::stopCluster(cl)
  if(permutation_rounds > 0){
    determine_fdr(output_file_name_cis = output_file_name_cis_ct_cond, permutation_rounds = permutation_rounds, do_smallest_per_gene = do_smallest_per_gene)
  }
}


perform_qtl_mapping <- function(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = 0, permute_in_covar_group=NULL, maf=0.1, cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('Baseline', 't24h', 't8w'), do_smallest_per_gene=T){
  # read the snps
  snps <- fread(snps_loc, sep = '\t', header = T, stringsAsFactors=FALSE)
  # read the positions of the gene and snp positions
  snpspos = fread(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  genepos = fread(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  # confinement of genes
  gene_confinement <- NULL
  # filter by confinement
  if(!is.null(confinement_file_name)){
    # read the table
    confinement <- read.table(confinement_file_name, sep = '\t', header = F, stringsAsFactors = F)
    # we check SNPs for sure
    snps_confined <- confinement$V1
    # and subset the SNPs
    snps <- snps[!is.na(snps$id) & snps$id %in% snps_confined, ]
    # check the number of columns
    if(ncol(confinement) > 1){
      # get the genes
      gene_confinement <- confinement$V2
    }
  }
  for(cell_type in cell_typers){
    for(condition in conditions){
      # get the specific expression data
      features_loc_ct_cond <- paste(features_loc_prepend, condition, '/',  cell_type, features_loc_append, sep = '') 
      output_file_name_cis_ct_cond <- paste(output_file_name_cis_prepend, condition, '/',  cell_type, output_file_name_cis_append, sep = '')
      covariates_file_loc_ct_cond <- paste(covariates_file_name_prepend, condition, '/',  cell_type, covariates_file_name_append, sep = '')
      # do the mapping
      run_qtl_mapping(features_loc_ct_cond, output_file_name_cis_ct_cond, covariates_file_loc_ct_cond, snps, snpspos, genepos, maf, permutation_rounds = permutation_rounds , permute_in_covar_group = permute_in_covar_group, do_smallest_per_gene = do_smallest_per_gene, gene_confinement=gene_confinement)
    }
  }
}

setup_default_settings <- function(){
  # set some defaults here
  default_settings_list <- list(
    maf = 0.1,
    permutation_rounds = 20,
    features_loc_append = '_expression.tsv',
    output_file_name_cis_append = '.cis.tsv',
    permute_in_covar_group = 'chem',
    cell_types=paste(c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), collapse = ',')
  )
  # return these
  return(default_settings_list)
}

setup_settings <- function(configuration_table){
  # first load the defaults
  configuration_list <- setup_default_settings()
  # start overwriting the defaults
  for(i in 1:nrow(configuration_table)){
    # add these rows, which are key-value pairs, to the config
    configuration_list[[configuration_table[i, 1]]] <- configuration_table[i, 2]
  }
  # return the configuration list
  return(configuration_list)
}

# do CLI method, when used as script
cli <- F
if(cli){
  # read the command line arguments
  args <- commandArgs(trailingOnly=TRUE)
  # the config file is a parameter
  config_loc <- args[1]
  # read the file
  config_table <- read.table(config_loc, header = F, sep = '\t')
  # set up the configuration
  settings_list <- setup_settings(config_table)
  
  # check for any confinements
  confinement_file_name <- NULL
  if('confinement_file_name' %in% names(settings_list)){
    confinement_file_name <- settings_list[['confinement_file_name']]
  }
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  snps_loc <- settings_list[['snps_loc']]
  #snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
  snps_location_file_name <- settings_list[['snps_location_file_name']]
  #snps_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name <- settings_list[['gene_location_file_name']]
  #gene_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend <- settings_list[['features_loc_prepend']]
  #features_loc_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_all_lowerres_20210301_metaqtl/'
  output_file_name <- settings_list[['output_file_name']]
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_harst_100k/'
  covariates_file_name_prepend <- settings_list[['covariates_file_name_prepend']]
  #covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_all_lowerres_20210301/'
  
  features_loc_append <- settings_list[['features_loc_append']]
  #features_loc_append<-'_expression.tsv'
  output_file_name_cis_append <- settings_list[['output_file_name_cis_append']]
  #output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append <- settings_list[['covariates_file_name_append']]
  #covariates_file_name_append<-'_metadata.chem.tsv'
  
  maf <- settings_list[['maf']]
  #maf <- 0.1
  permutation_rounds <- settings_list[['permutation_rounds']]
  #permutation_rounds <- 20
  
  
  cell_typers <- strsplit(settings_list[['cell_types']], ',')[[1]]
  #cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- strsplit(settings_list[['cell_types']], ',')[[1]]
  #conditions <- c('Baseline', 't24h', 't8w')
  
  permute_in_covar_group <- settings_list[['permute_in_covar_group']]
  #permute_in_covar_group <- 'chem'
  
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}

do_all_stemi_meta <- F
if(do_all_stemi_meta){
  confinement_file_name <- NULL
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_all_lowerres_20210301_metaqtl/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_harst_100k/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_eqtlgenlead/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_all_lowerres_20210301/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- c('Baseline', 't24h', 't8w')
  
  permute_in_covar_group <- 'chem'
  
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}

# UT+HC
do_all_ut_stemi <- F
if(do_all_ut_stemi){
  confinement_file_name <- NULL
  confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_harst_100k/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_eqtlgenlead/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem_status.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w', 'UT')
  
  permute_in_covar_group <- 'chem'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}


all_env <- NULL


do_all_ut_stemi_eqtlgen <- F
if(do_all_ut_stemi_eqtlgen){
  confinement_file_name <- NULL
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  confinement_file_name <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_harst_100k/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_eqtlgenlead/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_eqtlgenlead/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem_status.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w', 'UT')
  
  permute_in_covar_group <- 'chem'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}

do_all_ut_stemi_unconfined <- T
if(do_all_ut_stemi_unconfined){
  confinement_file_name <- NULL
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_unconfined_100k/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem_status.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w', 'UT')
  
  permute_in_covar_group <- 'chem'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}

