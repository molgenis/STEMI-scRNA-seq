# load libraries
library(Seurat)
library(Matrix)
library(ggplot2)

plot_qtl <- function(mean_expression_matrix, genotypes, gene, snp){
  # grab the common participants
  common_participants <- intersect(colnames(mean_expression_matrix), colnames(genotypes))
  # grab the expression
  expression <- as.vector(unlist(mean_expression_matrix[gene, common_participants]))
  # grab the SNPs
  genotype <- as.vector(unlist(genotypes[snp, common_participants]))
  # create a dataframe
  plot_frame <- data.frame(genotype=genotype, expression=expression)
  # make the plot
  p <- ggplot(data=plot_frame, mapping=aes(x=genotype, y=expression, fill=genotype)) + geom_boxplot()
  return(p)
}

mash_matrixeqtl_output <- function(base_output_path, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('Baseline', 't24h','t8w')){
  
}

get_conditions_mash <- function(eqtl_output_loc, conditions = c('Baseline', 't24h','t8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  betas = NULL
  ses = NULL
  for(cell_type in cell_types){
    # go through the conditions
    for(condition in conditions){
      # there might be some tables not present, so let's set to a default
      table <- NULL
      # read the table
      table_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/', 'eQTLsFDR-ProbeLevel.txt.gz',  sep = '')
      # show progress
      print(paste('reading', table_loc))
      try({
        table <- read.table(table_loc, sep = '\t', header = T)
      })
      if(is.null(table)){
        print(paste('table not present for:', table_loc))
      } else{
        # grab the variables
        table_b1 <- as.numeric(str_match(table$Beta..SE., "(.+) \\x28.+\\x29;")[,2])
        table_se1 <- as.numeric(str_match(table$Beta..SE., ".+ \\x28(.+)\\x29;")[,2])
        table_b2 <- as.numeric(str_match(table$Beta..SE., ";(.+) \\x28.+\\x29")[,2])
        table_se2 <- as.numeric(str_match(table$Beta..SE., ";.+ \\x28(.+)\\x29")[,2])
        # add the SNP>probe as rownames
        table_b1 <- data.frame(table_b1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        table_se1 <- data.frame(table_se1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        table_b2 <- data.frame(table_b2, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        table_se2 <- data.frame(table_se2, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        # add the colnames
        colnames(table_b1) <- c(paste(condition, 'v2', sep = '.'))
        colnames(table_se1) <- c(paste(condition, 'v2', sep = '.'))
        colnames(table_b2) <- c(paste(condition, 'v3', sep = '.'))
        colnames(table_se2) <- c(paste(condition, 'v3', sep = '.'))
        # banana
        if(is.null(ses)){
          ses <- merge(table_se1, table_se2, by="row.names",all.x=TRUE)
          rownames(ses) <- ses$Row.names
          ses$Row.names <- NULL
        } else{
          ses <- merge(ses, table_se1, by="row.names",all.x=TRUE)
          rownames(ses) <- ses$Row.names
          ses$Row.names <- NULL
          ses <- merge(ses, table_se2, by="row.names",all.x=TRUE)
          rownames(ses) <- ses$Row.names
          ses$Row.names <- NULL
        }
        if(is.null(betas)){
          betas <- merge(table_b1, table_b2, by="row.names",all.x=TRUE)
          rownames(betas) <- betas$Row.names
          betas$Row.names <- NULL
        } else{
          betas <- merge(betas, table_b1, by="row.names",all.x=TRUE)
          rownames(betas) <- betas$Row.names
          betas$Row.names <- NULL
          betas <- merge(betas, table_b2, by="row.names",all.x=TRUE)
          rownames(betas) <- betas$Row.names
          betas$Row.names <- NULL
        }
      }
    }
  }
  # in case reading failed, or we have only one column we won't continue
  data  <- NULL
  if(!is.null(ses) & !is.null(betas) & ncol(ses) > 1 & ncol(betas) > 1){
    # put it in a list
    ses <- ses[complete.cases(ses), ]
    betas <- betas[complete.cases(betas), ]
    data = mash_set_data(as.matrix(betas), as.matrix(ses))
  }
  # return the mash data set
  return(data)
}

get_conditions_mash_meta_z <- function(eqtl_output_loc, conditions = c('Baseline', 't24h','t8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  betas = NULL
  ses = NULL
  # go through the cell types
  for(cell_type in cell_types){
    # go through the conditions
    for(condition in conditions){
      # there might be some tables not present, so let's set to a default
      table <- NULL
      # read the table
      table_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/', 'eQTLsFDR-ProbeLevel.txt.gz',  sep = '')
      # show progress
      print(paste('reading', table_loc))
      try({
        table <- read.table(table_loc, sep = '\t', header = T)
      })
      if(is.null(table)){
        print(paste('table not present for:', table_loc))
      } else{
        # grab the variables
        table_b1 <- table$OverallZScore
        # add the SNP>probe as rownames
        table_b1 <- data.frame(table_b1, row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        # add the colnames
        colnames(table_b1) <- c(condition)
        # create the SEs
        table_se1 <- data.frame(c(rep(1, nrow(table_b1))), row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        colnames(table_se1) <- c(condition)
        # banana
        if(is.null(ses)){
          ses <- table_se1
        } else{
          ses <- merge(ses, table_se1, by="row.names",all.x=TRUE)
          rownames(ses) <- ses$Row.names
          ses$Row.names <- NULL
        }
        if(is.null(betas)){
          betas <- table_b1
        } else{
          betas <- merge(betas, table_b1, by="row.names",all.x=TRUE)
          rownames(betas) <- betas$Row.names
          betas$Row.names <- NULL
        }
      }
    }
  }
  # in case reading failed, or we have only one column we won't continue
  data  <- NULL
  if(!is.null(ses) & !is.null(betas) & ncol(ses) > 1 & ncol(betas) > 1){
    # put it in a list
    ses <- ses[complete.cases(ses), ]
    betas <- betas[complete.cases(betas), ]
    data = mash_set_data(as.matrix(betas), as.matrix(ses))
  }
  # return the mash data set
  return(data)
}


get_conditions_cell_types_mash <- function(eqtl_output_loc, conditions = c('Baseline', 't24h','t8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_z=F){
  data <- NULL
  if(use_z){
    data <- get_conditions_mash_meta_z(eqtl_output_loc, conditions = conditions, cell_types)
  } else{
    data <- get_conditions_mash(eqtl_output_loc, conditions = conditions, cell_types)
  }
  lfsr <- NULL
  # result might be NULL if there is not enough data to compare
  if(!is.null(data)){
    # run MASH
    m.1by1 = mash_1by1(data)
    strong <- get_significant_results(m.1by1, 0.05)
    U.pca = cov_pca(data,ncol(data$Bhat), subset = strong) # 6PCS due to number of conditions
    U.ed = cov_ed(data, U.pca)
    U.c = cov_canonical(data)
    m = mash(data, c(U.c,U.ed))
    lfsr <- m$result$lfsr
  }
  return(lsfr)
}


# load object
#cardio.stemi <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.stemi.20210611.combatcorrected.rds')
# genotype location
snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
# load genotype data
snps <- fread(snps_loc, header=T)
# locations of features
features_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/'


