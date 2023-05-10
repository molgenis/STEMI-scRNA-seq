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

determine_fdr <- function(output_file_name_cis, permutation_rounds, do_smallest_per_gene=T, use_datatable=T, verbose_number=NA){
  # we will add the permutation rounds together
  permutation_table <- NULL
  # we'll store the SNP and gene info separately, that's easier for the rowmeans
  snp_gene <- NULL
  for(i in 1:permutation_rounds){
    if(use_datatable){
      # read the file
      permutation_round <- fread(paste(output_file_name_cis, '.permuted.', i, sep = ''), sep = '\t', header = T, check.names = T) # check.names T to keep compatibility between data.frame and data.table
      # remove columns we do not need
      permutation_round[, beta:=NULL]
      permutation_round[, t.stat:=NULL]
      permutation_round[, FDR:=NULL]
      # rename the p.value column
      permutation_round[[paste('perm', i, sep='')]] <- permutation_round$p.value
      permutation_round[, p.value:=NULL]
      # set the key
      setkey(permutation_round, SNP, gene)
      # check if the first round
      if(is.null(permutation_table)){
        permutation_table <- permutation_round
      }
      else{
        permutation_table <- merge(permutation_table, permutation_round, by=c('SNP', 'gene'), all=T)
      }
    }
    else{
      # reset the order
      permutation_round <- read.table(paste(output_file_name_cis, '.permuted.', i, sep = ''), sep = '\t', header = T, stringsAsFactors = F)
      rownames(permutation_round) <- paste(permutation_round[['SNP']], permutation_round[['gene']], sep = '_')
      # paste onto the rest of the output
      if(is.null(permutation_table)){
        #permutation_table <- data.frame(perm1 = permutation_round[['p-value']])
        permutation_table <- permutation_round[, c('p.value'), drop = F]
        permutation_table[['perm1']] <- permutation_round[['p.value']]
        permutation_round[['p.value']] <- NULL
        # we'll need this one later
        snp_gene <- permutation_round[, c('SNP', 'gene')]
      }
      else{
        # grab the same snp-gene pairs
        permutation_table[[paste('perm', i, sep = '')]] <- permutation_round[rownames(permutation_table), c('p.value'), drop = F]
      }
    }
  }
  permuted_table <- NULL
  print('creating mean permuted p per snp-probe')
  if(use_datatable){
    # fetch the SNP and gene columns
    snps <- permutation_table[['SNP']]
    genes <- permutation_table[['gene']]
    # remove columns that interfere with the rownmeans
    permutation_table[, SNP:=NULL]
    permutation_table[, gene:=NULL]
    # calculate the average p-value
    mean_p_permuted=rowMeans(permutation_table)
    # create data.table
    permuted_table <- data.table(SNP=snps, gene=genes, perm=mean_p_permuted)
    # order by significance
    print('sorting by mean permuted p')
    setorder(permuted_table, perm)
  }
  else{
    # calculate the average p-value
    mean_p_permuted=rowMeans(permutation_table)
    # add back the snp and gene info
    permuted_table <- snp_gene
    permuted_table[['perm']] <- mean_p_permuted
    # order by significance
    print('sorting by mean permuted p')
    permuted_table <- permuted_table[order(permuted_table[['perm']]), ]
  }
  # filter to the best value per gene, if requested
  if(do_smallest_per_gene){
    print('selecting lowest permutated p per gene')
    # we ordered the table, so the first entry found for a gene, should always be the most significant one
    permuted_table <- permuted_table[match(unique(permuted_table[['gene']]), permuted_table[['gene']]), ]
  }
  # now read the non-permuted file
  output_file <- fread(paste(output_file_name_cis), sep = '\t', header = T, check.names = T)
  setorder(output_file, p.value)
  #output_file <- read.table(paste(output_file_name_cis), sep = '\t', header = T, stringsAsFactors = F)
  # go through each row
  print('determening permuted FDR for each SNP-probe combination')
  #output_file[['permuted_fdr']] <- apply(output_file, 1, function(x){
  #  # grab the p-value
  #  p_real <- as.numeric(x[['p.value']])
  #  # skip the ones that can't possibly be significant
  #  if(p_real == 1){
  #    return(1)
  #  }
  #  else{
  #    # get the number of permuted p values
  #    nr_permuted_p <- nrow(permuted_table)
  #    # check how many of these p values are smaller than the actual p
  #    #nr_permuted_p_smaller <- nrow(permuted_table[permuted_table[['perm']] <= p_real, ])
  #    # get the first index where a permuted variable is as good or better than the actual p
  #    first_better_perm_p <- get_first_index_at_cutoff(permuted_table[['perm']], p_real)
  #    # substract one, because the indexing starts at 1, but that is the first time the evaluation failed
  #    nr_permuted_p_smaller <- first_better_perm_p - 1
  #    # calculate the fraction of permuted Ps, that are smaller than your actual p
  #    perm_fdr <- nr_permuted_p_smaller / nr_permuted_p
  #    # that fraction is the chance that your p could come from the permuted distribution
  #    return(perm_fdr)
  #  }
  #})
  # create vector to hold the calculated FDRs
  permuted_fdr <- vector(mode='numeric', length=nrow(output_file))
  # grab only the permuted p-values
  permuted_ps <- permuted_table$perm
  # get the real ps
  real_ps <- output_file$p.value
  # because our actual variables are also sorted, we can keep track of where we were in the permuted set as well
  index_permuted <- 1
  # check each variable
  for(p_real_index in 1:length(real_ps)){
    # set when to stop
    found_higher_perm <- F
    # grab the actual p value
    p_real <- real_ps[p_real_index]
    # check each permuted variable whether or not it is better
    while(found_higher_perm == F & index_permuted <= length(permuted_ps)){
      # grab the permuted p
      permuted_p <- permuted_ps[index_permuted]
      # check if the same
      if(permuted_p == p_real){
        # go until you see a new P value
        if(index_permuted != length(permuted_ps) & permuted_p == permuted_ps[index_permuted+1]){
          # we can just skip to the next one
          index_permuted <- index_permuted+1
        }
        # if the next is a new value, this is the number that was better or equal
        nr_permuted_equal_or_better <- index_permuted
        # stop search here
        found_higher_perm <- T
        # add result to the fdrs
        permuted_fdr[p_real_index] <- nr_permuted_equal_or_better/length(permuted_ps)
        # stop searching and move on to the next p
        found_higher_perm <- T
      }
      # if the permuted p is worse than the original p
      if(permuted_p > p_real){
        # the last p value we saw was the last one that was better
        nr_permuted_equal_or_better <- index_permuted - 1
        # stop search here
        found_higher_perm <- T
        # add result to the fdrs
        permuted_fdr[p_real_index] <- nr_permuted_equal_or_better/length(permuted_ps)
        # stop searching and move on to the next p
        found_higher_perm <- T
      }
      # otherwise we will move on to check the next permuted p
      else{
        index_permuted <- index_permuted+1
      }
    }
    if(found_higher_perm == F){
      permuted_fdr[p_real_index] <- 1
    }
    # both vectors are sorted by their size, as such we can keep continueing to search from where we left in the permuted vector, as the next real p will be of equal size or bigger
  }
  # add the FDRs we calculated
  output_file$permuted_fdr <- permuted_fdr
  # also add bonferoni
  output_file$bonferroni <- p.adjust(output_file$p.value, method = 'bonferroni')
  # get the file without the extention
  output_file_name_cis_fdred <- sub('\\..[^\\.]*$', '', output_file_name_cis)
  # add new extention
  output_file_name_cis_fdred <- paste(output_file_name_cis_fdred, '.fdr.tsv', sep = '')
  # write the file
  print('writing results with permuted FDR')
  write.table(output_file, output_file_name_cis_fdred, sep = '\t', row.names = F, col.names = T, quote = F)
}

determine_fdr_emp <- function(output_file_name_cis, permutation_rounds, verbose_number=NA){
  # read the unpermuted file
  output_file <- fread(paste(output_file_name_cis), sep = '\t', header = T, check.names = T)
  # sort the p values
  setorder(output_file, p.value)
  # get the unique genes, this will be the order we will use when the order matters
  genes <- unique(output_file[['gene']])
  # grab the top effect per gene, because it is ordered, the first match is the one with the lowest p
  ps_real <- as.vector(unlist(output_file[match(genes, output_file[['gene']]), 'p.value']))
  # now reserve space for the p values
  ps_permuted <- vector(mode = 'numeric', length = (length(genes) * permutation_rounds))
  for(i in 1:permutation_rounds){
    # read the file
    permutation_round <- fread(paste(output_file_name_cis, '.permuted.', i, sep = ''), sep = '\t', header = T, check.names = T) # check.names T to keep compatibility between data.frame and data.table
    # sort by p value
    setorder(permutation_round, p.value)
    # grab the top effect per gene, because it is ordered, the first match is the one with the lowest p
    ps_permuted_round <- as.vector(unlist(permutation_round[match(genes, permutation_round[['gene']]), 'p.value']))
    # determine where to fill the vector
    end_pos <- i*length(genes)
    start_pos <- ((i-1)*length(genes)) + 1
    # add the p values
    ps_permuted[start_pos:end_pos] <- ps_permuted_round
  }
  # sort the permuted p values
  ps_permuted <- ps_permuted[order(ps_permuted)]
  # we will have an FDR for each gene, so let's reserve space once again
  gene_fdrs <- vector(mode = 'numeric', length = length(ps_real))
  # check each variable
  for(p_real_index in 1:length(ps_real)){
    # check the proportion of P values that is better than in the permutations
    nr_ps_permuted_better <- sum(ps_permuted <= ps_real[p_real_index])
    # check how many were better than the actual ps, correcting for the number of permutations
    fdr_gene <- nr_ps_permuted_better/length(genes)/permutation_rounds
    # put this in our gene-based FDR vector
    gene_fdrs[p_real_index] <- fdr_gene
    # report on progress if that makes one happy
    if(!is.na(verbose_number) & p_real_index %% verbose_number == 0){
      print(paste('checked', p_real_index, 'genes'))
    }
  }
  # put the gene and FDR in a datatable, remember we did this in the same order
  fdr_table <- data.table(gene=genes, fdr_gene=gene_fdrs)
  setkey(fdr_table, gene)
  # merge with the output file we already have, datatables merge very efficiently
  #setkey(output_file, cols=c('gene', 'SNP'))
  output_file$FDR_gene <- as.vector(unlist(fdr_table[match(output_file[['gene']], fdr_table[['gene']]), 'fdr_gene']))
  # also add bonferoni
  output_file$bonferroni <- p.adjust(output_file$p.value, method = 'bonferroni')
  # get the file without the extention
  output_file_name_cis_fdred <- sub('\\..[^\\.]*$', '', output_file_name_cis)
  # add new extention
  output_file_name_cis_fdred <- paste(output_file_name_cis_fdred, '.fdr.tsv', sep = '')
  # write the file
  print('writing results with permuted gene FDR')
  write.table(output_file, output_file_name_cis_fdred, sep = '\t', row.names = F, col.names = T, quote = F)
}

run_qtl_mapping <- function(features_loc_ct_cond, output_file_name_cis_ct_cond, covariates_file_loc_ct_cond, snps, snpspos, genepos, maf=0.1, permutation_rounds = 0, permute_in_covar_group=NULL, do_smallest_per_gene=T, gene_confinement=NULL, snp_confinement=NULL, snp_gene_confinement=NULL){
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
  if(!is.null(snp_gene_confinement)){
    snps <- snps[!is.na(snps$id) & snps$id %in% snp_gene_confinement[[1]], ]
    expressions_ct_cond <- expressions_ct_cond[!is.na(expressions_ct_cond$id) & expressions_ct_cond$id %in% snp_gene_confinement[[2]], ]
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
  # calculate the maf for the SNPs that are left
  maf_calculated <- as.vector(unlist(rowSums(snps_ct_cond[, 2:ncol(snps_ct_cond)])/(ncol(snps_ct_cond)-1)/2))
  # grab the SNPs as well
  snps_to_be_tested <- snps_ct_cond$id
  # and the number of participants we had
  nr_participants <- ncol(snps_ct_cond) -1 # first column is the ID column
  # put into data table
  snps_metadata <- data.table(SNP=snps_to_be_tested, maf=maf_calculated, n=nr_participants)
  
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
  if(!is.null(snp_gene_confinement)){
    # because MatrixEQTL doesn't have snp+gene confinements, we need to do filtering after the fact
    result <- read.table(output_file_name_cis_ct_cond, header = T, sep = '\t', stringsAsFactors = F, check.names = F)
    # pasting SNP and gene together, allows us to filter
    result$snp_probe <- paste(result$SNP, result$gene, sep = '_')
    result <- result[result$snp_probe %in% paste(snp_gene_confinement[[1]], snp_gene_confinement[[2]], sep = '_'), ]
    # remove the column used for this filtering
    result$snp_probe <- NULL
    # recalculate the fdr that MatrixEQTL does (n is different after filtering)
    result$FDR <- p.adjust(result[['p-value']], method = 'fdr')
    # rewrite the result
    write.table(result, output_file_name_cis_ct_cond, sep = '\t', row.names = F, col.names = T) 
  }
  # add the SNP metadata back we created before the mapping
  result <- fread(output_file_name_cis_ct_cond, check.names = T)
  # key the snp and genek
  setkey(result)
  # merge the data
  result <- merge(result, snps_metadata, by='SNP')
  # calculate the r
  result$r <- result$t.stat / sqrt(result$n - 2 + result$t.stat ** 2)
  # and write our result with the metadata added
  write.table(result, output_file_name_cis_ct_cond, sep = '\t', row.names = F, col.names = T) 
  
  # convert data tables to data frames
  expressions_ct_cond <- data.frame(expressions_ct_cond)
  covariates_ct_cond <- data.frame(covariates_ct_cond)
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
      print(paste('doing grouped permutations with', permute_in_covar_group, sep = ' '))
      # start with the unpermuted file
      expressions_ct_cond_permuted <- expressions_ct_cond
      # covariates as well
      covariates_ct_cond_permuted <- covariates_ct_cond
      # get the groups to permute in
      vars_in_group <- as.vector(unlist(covariates_ct_cond[covariates_ct_cond[['id']] == permute_in_covar_group, 2:ncol(covariates_ct_cond)])) # skipping column one, because that is the 'row name'
      # go through each permutable category
      for(group in unique(vars_in_group)){
        print(paste('permuting for group', group, sep = ' '))
        # get the samples with this group
        participant_locations <- which(vars_in_group == group) + 1 # adding one to shift right, to account for skipping the 'row name' column
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
      print('doing ungrouped permutations')
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
    if(!is.null(snp_gene_confinement)){
      # because MatrixEQTL doesn't have snp+gene confinements, we need to do filtering after the fact
      perm_result <- read.table(output_file_name_cis_ct_cond_permuted, header = T, sep = '\t', stringsAsFactors = F, check.names = F)
      # pasting SNP and gene together, allows us to filter
      perm_result$snp_probe <- paste(perm_result$SNP, perm_result$gene, sep = '_')
      perm_result <- perm_result[perm_result$snp_probe %in% paste(snp_gene_confinement[[1]], snp_gene_confinement[[2]], sep = '_'), ]
      # remove the column used for this filtering
      perm_result$snp_probe <- NULL
      # recalculate the fdr that MatrixEQTL does (n is different after filtering)
      perm_result$FDR <- p.adjust(perm_result[['p-value']], method = 'fdr')
      # rewrite the result
      write.table(perm_result, output_file_name_cis_ct_cond_permuted, sep = '\t', row.names = F, col.names = T)  
    }
  }
  #parallel::stopCluster(cl)
  if(permutation_rounds > 0){
    determine_fdr_emp(output_file_name_cis = output_file_name_cis_ct_cond, permutation_rounds = permutation_rounds)
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
  # confinement of snp_gene combinations
  snp_gene_confinement <- NULL
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
      # and set the combination of snp-gene confinement
      snp_gene_confinement <- confinement
    }
  }
  for(cell_type in cell_typers){
    for(condition in conditions){
      # get the specific expression data
      features_loc_ct_cond <- paste(features_loc_prepend, condition, '/',  cell_type, features_loc_append, sep = '') 
      output_file_name_cis_ct_cond <- paste(output_file_name_cis_prepend, condition, '/',  cell_type, output_file_name_cis_append, sep = '')
      covariates_file_loc_ct_cond <- paste(covariates_file_name_prepend, condition, '/',  cell_type, covariates_file_name_append, sep = '')
      # do the mapping
      run_qtl_mapping(features_loc_ct_cond, output_file_name_cis_ct_cond, covariates_file_loc_ct_cond, snps, snpspos, genepos, maf, permutation_rounds = permutation_rounds , permute_in_covar_group = permute_in_covar_group, do_smallest_per_gene = do_smallest_per_gene, gene_confinement=gene_confinement, snp_gene_confinement = snp_gene_confinement)
    }
  }
}


re_annotate_qlt_output <- function(qtl_output_loc, ref_table_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w'), sort_columns=c('FDR_gene', 'p.value')){
  # read the reference table
  ref_table <- fread(ref_table_loc)
  # set the same column as in the eQTL mapping output
  #ref_table$SNP <- ref_table$ID
  #ref_table$ID <- NULL
  # set as key for fast mergine
  setkey(ref_table, SNP)
  # for each cell type and condition
  for(cell_type in cell_types){
    for(condition in conditions){
      try({
        # read the original table
        eqtl_table <- fread(paste(qtl_output_loc, condition, '/', cell_type, '.cis.fdr.tsv', sep = ''))
        # also set key for fast merging
        setkey(eqtl_table, SNP)
        # merge together
        eqtl_table <- merge(eqtl_table, ref_table, by = 'SNP', all = F)
        output_loc_annotated <- paste(qtl_output_loc, condition, '/', cell_type, '.cis.fdr.gt.tsv', sep = '')
        print(head(eqtl_table))
        # re-sort on gene FDR
        #setorder(eqtl_table, cols=c(sort_columns))
        # write result
        print(paste('writing to', output_loc_annotated))
        write.table(eqtl_table, output_loc_annotated, sep = '\t', row.names = F, col.names = T, quote = F)
      })
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

do_all_stemi_meta <- T
if(do_all_stemi_meta){
  confinement_file_name <- NULL
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_all_lowerres_20210301_metaqtl/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_harst_100k/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_eqtlgenlead/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_all_lowerres_20210301/'
  
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
do_all_ut_stemi <- T
if(do_all_ut_stemi){
  confinement_file_name <- NULL
  confinement_file_name <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_harst_100k/'
  #output_file_name_cis_prepend<-'/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301_eqtlgenlead/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem_status.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK', 'bulk')
  conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w')
  
  permute_in_covar_group <- 'chem_V3'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}


all_env <- NULL


do_all_ut_stemi_eqtlgen <- T
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
  
  cell_typers=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w')
  
  permute_in_covar_group <- 'chem_V3'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
}

do_all_ut_stemi_unconfined <- T
if(do_all_ut_stemi_unconfined){
  confinement_file_name <- NULL
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/harst_2017_SNPs.txt'
  #confinement_file_name <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/confinements/eqtl_v1013_lead_snp_gene.txt'
  
  snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric.tsv'
  #snps_loc<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/hc_and_stemi_numeric_firstfour.tsv'
  snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
  #snps_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos_firstfour.tsv'
  gene_location_file_name<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'
  
  features_loc_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_and_1mut_lowerres_20210629_metaqtl/'
  output_file_name_cis_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_and_1mut_lowerres_20210629_unconfined_100k/'
  covariates_file_name_prepend<-'/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_and_1mut_lowerres_20210629/'
  
  features_loc_append<-'_expression.tsv'
  output_file_name_cis_append<-'.cis.tsv'
  covariates_file_name_append<-'_metadata.chem_status.tsv'
  
  maf <- 0.1
  permutation_rounds <- 20
  
  #cell_typers=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
  #conditions <- c('UT_Baseline', 'UT_t24h', 'UT_t8w')
  
  permute_in_covar_group <- 'chem_V3'
  #snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append
  #perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
  
  cell_typers=c('B')
  conditions <- c('UT_t24h', 'UT_t8w')
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
  cell_typers=c('CD8T')
  conditions <- c('UT_t8w')
  perform_qtl_mapping(snps_loc, snps_location_file_name, gene_location_file_name, features_loc_prepend, output_file_name_cis_prepend, covariates_file_name_prepend, features_loc_append, output_file_name_cis_append, covariates_file_name_append, confinement_file_name, permutation_rounds = permutation_rounds, permute_in_covar_group=permute_in_covar_group, cell_typers=cell_typers, conditions=conditions)
  
}

