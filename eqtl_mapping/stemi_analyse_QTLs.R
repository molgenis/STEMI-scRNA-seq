# load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(UpSetR)

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
        colnames(table_b1) <- c(paste(condition, cell_type, 'v2', sep = '.'))
        colnames(table_se1) <- c(paste(condition, cell_type, 'v2', sep = '.'))
        colnames(table_b2) <- c(paste(condition, cell_type, 'v3', sep = '.'))
        colnames(table_se2) <- c(paste(condition, cell_type, 'v3', sep = '.'))
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
        colnames(table_b1) <- c(paste(cell_type, '.', condition, sep = ''))
        # create the SEs
        table_se1 <- data.frame(c(rep(1, nrow(table_b1))), row.names = paste(table$SNPChr,table$SNPChrPos,table$ProbeName,sep="_"))
        colnames(table_se1) <- c(paste(cell_type, '.', condition, sep = ''))
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


get_conditions_cell_types_mash <- function(eqtl_output_loc, conditions = c('Baseline', 't24h','t8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_z=F, npcs=5, full=F){
  data <- NULL
  if(use_z){
    data <- get_conditions_mash_meta_z(eqtl_output_loc, conditions = conditions, cell_types)
  } else{
    data <- get_conditions_mash(eqtl_output_loc, conditions = conditions, cell_types)
  }
  m <- NULL
  # result might be NULL if there is not enough data to compare
  if(!is.null(data)){
    # run MASH
    print('running mash')
    m.1by1 = mash_1by1(data)
    print('getting significant results')
    strong <- get_significant_results(m.1by1, 0.05)
    print('getting PCA')
    U.pca = cov_pca(data, npcs, subset = strong) # 5PCS due to number of conditions
    print('cov ed')
    U.ed = cov_ed(data, U.pca)
    print('U.c')
    U.c = cov_canonical(data)
    print('mashing')
    m = mash(data, c(U.c,U.ed))
  }
  if(full){
    return(m)
  }
  else{
    lfsr <- m$result$lfsr
    return(lfsr)
  }
}

check_exclusivity <- function(list_of_vectors, exclusivity_list){
  # extract the entry of the list to check for exclusivity
  exclusivity_vector <- list_of_vectors[[exclusivity_list]]
  # get the names of the other entries
  other_vectors_names <- setdiff(names(list_of_vectors), exclusivity_list)
  # grab the entries of the other vectors, and flatten them into one vector
  other_vector_entries <- unlist(list_of_vectors[other_vectors_names])
  # the extract what is only found in our vector to check for exclusivity
  exclusives <- setdiff(exclusivity_vector, other_vector_entries)
  # check if there are any, and return the correct answer
  if(length(exclusives) > 0){
    return(T)
  }
  else{
    return(F)
  }
}



get_matrixeqtl_egenes_per_cell_type <- function(output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene'){
  # 
  egenes_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # build the full eQTL path
    eqtl_path <- paste(output_loc, cell_type, '.cis.fdr.tsv', sep = '')
    # read the output file
    eqtl_output <- read.table(eqtl_path, sep = '\t', header = T, stringsAsFactors = F)
    # subset to the egenes
    eqtl_output_psig <- eqtl_output[eqtl_output[[pval_column]] < pval_cutoff, ]
    # if requested, check the gene FDR
    if(only_gene_fdr){
      eqtl_output_psig <- eqtl_output_psig[eqtl_output_psig[[gene_fdr_column]] < pval_cutoff, ]
    }
    # add the significant results to the list of this specific cell type
    if(nrow(eqtl_output_psig) > 0){
      egenes_per_ct[[cell_type]] <- unique(eqtl_output_psig[['gene']])
    }
  }
  return(egenes_per_ct)
}


get_matrixeqtl_egenes_per_condition <- function(output_loc, append, conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene'){
  # 
  egenes_per_condition <- list()
  # check each cell type
  for(condition in conditions){
    # build the full eQTL path
    eqtl_path <- paste(output_loc, condition, '/', append, sep = '')
    # read the output file
    eqtl_output <- read.table(eqtl_path, sep = '\t', header = T, stringsAsFactors = F)
    # subset to the egenes
    eqtl_output_psig <- eqtl_output[eqtl_output[[pval_column]] < pval_cutoff, ]
    # if requested, check the gene FDR
    if(only_gene_fdr){
      eqtl_output_psig <- eqtl_output_psig[eqtl_output_psig[[gene_fdr_column]] < pval_cutoff, ]
    }
    # add the significant results to the list of this specific cell type
    if(nrow(eqtl_output_psig) > 0){
      egenes_per_condition[[condition]] <- unique(eqtl_output_psig[['gene']])
    }
  }
  return(egenes_per_condition)
}


plot_egene_sharing_per_celltype <- function(output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene', use_label_dict=T, use_color_dict=T){
  # get the 
  egenes_per_ct <- get_matrixeqtl_egenes_per_cell_type(output_loc, cell_types=cell_types, pval_column=pval_column, pval_cutoff=pval_cutoff, only_gene_fdr=only_gene_fdr, gene_fdr_column=gene_fdr_column)
  # change the labels to something more pretty
  if(use_label_dict){
    names(egenes_per_ct) <- label_dict()[names(egenes_per_ct)]
  }
  # if colouring this list is filled
  queries <- NULL
  # this parameter must be supplied due to a bug
  empty.intersections <-  NULL
  # default color
  sets.bar.color <- 'black'
  # upset has a bug with empty intersections
  dont_colour <- F
  # set coloured bars for the wholly unique cell types
  if(use_color_dict){
    queries <- list()
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    for(i in 1:length(names(egenes_per_ct))){
      # the name is the cell type
      cell_type <- names(egenes_per_ct)[i]
      # add for the singles in the intersection sizes
      ct_list <- list(
        query = intersects,
        params = list(cell_type),
        color = get_color_coding_dict()[[cell_type]],
        active = T)
      queries[[i]] <- ct_list
      # check for exclusivity, there is a bug with colouring of empty intersections
      is_exclusive <- check_exclusivity(egenes_per_ct, cell_type)
      if(!is_exclusive){
        empty.intersections <- "on"
        # we can't colour with empty intersections unfortunately
        dont_colour <- T
      }
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(cell_type), nr=c(length(egenes_per_ct[[cell_type]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
    # add the colors for the sets
    sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
  }
  # if there are empty queries, we sadly have to forego the colouring, as there is but in upsetr with these
  if(dont_colour){
    queries <- NULL
  }
  return(upset(fromList(egenes_per_ct), order.by = 'freq', nsets = length(names(egenes_per_ct)), queries = queries, sets.bar.color=sets.bar.color, empty.intersections = empty.intersections))
}


plot_egene_sharing_per_condition <- function(output_loc, append, conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene', use_label_dict=T, use_color_dict=T){
  # get the 
  egenes_per_condition <- get_matrixeqtl_egenes_per_condition(output_loc, append, conditions=conditions, pval_column=pval_column, pval_cutoff=pval_cutoff, only_gene_fdr=only_gene_fdr, gene_fdr_column=gene_fdr_column)
  # change the labels to something more pretty
  if(use_label_dict){
    names(egenes_per_condition) <- label_dict()[names(egenes_per_condition)]
  }
  # if colouring this list is filled
  queries <- NULL
  # this parameter must be supplied due to a bug
  empty.intersections <-  NULL
  # default color
  sets.bar.color <- 'black'
  # set coloured bars for the wholly unique cell types
  if(use_color_dict){
    queries <- list()
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    for(i in 1:length(names(egenes_per_condition))){
      # the name is the cell type
      condition <- names(egenes_per_condition)[i]
      # add for the singles in the intersection sizes
      ct_list <- list(
        query = intersects,
        params = list(condition),
        color = get_color_coding_dict()[[condition]],
        active = T)
      queries[[i]] <- ct_list
      # check for exclusivity, there is a bug with colouring of empty intersections
      is_exclusive <- check_exclusivity(egenes_per_condition, condition)
      if(is_exclusive){
        empty.intersections <- "on"
      }
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(condition), nr=c(length(egenes_per_condition[[condition]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
    # add the colors for the sets
    sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
  }
  return(upset(fromList(egenes_per_condition), order.by = 'freq', nsets = length(names(egenes_per_condition)), queries = queries, sets.bar.color=sets.bar.color, empty.intersections = empty.intersections))
}


plot_egene_sharing_per_condition_allct <- function(output_loc, append, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene', use_label_dict=T, use_color_dict=T){
  # make plot for all cell types
  plot_per_cell_type <- list()
  # go through the cell types
  for(cell_type in cell_types){
    # set up the append that we will supply
    ct_append <- paste(cell_type, append, sep = '')
    plot_ct <- plot_egene_sharing_per_condition(output_loc = output_loc, append = ct_append, conditions = conditions, pval_column = pval_column, pval_cutoff = pval_cutoff, only_gene_fdr = only_gene_fdr, gene_fdr_column = gene_fdr_column, use_label_dict = use_label_dict, use_color_dict = use_color_dict)
    # put the plot in the list
    plot_per_cell_type[[cell_type]] <- plot_ct
  }
  return(plot_per_cell_type)
}

plot_egene_sharing_per_celltype_allcond <- function(output_loc, conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), pval_column='p.value', pval_cutoff=0.05, only_gene_fdr=T, gene_fdr_column='FDR_gene', use_label_dict=T, use_color_dict=T){
  # make the plot for each condition
  plot_per_condition <- list()
  # go through the conditions
  for(condition in conditions){
    # set up the path for this specific condition
    output_loc_condition <- paste(output_loc, condition, '/', sep = '')
    # make the plots
    plot_condition <- plot_egene_sharing_per_celltype(output_loc = output_loc_condition, cell_types = cell_types, pval_column = pval_column, pval_cutoff = pval_cutoff, only_gene_fdr = only_gene_fdr, gene_fdr_column = gene_fdr_column, use_label_dict = use_label_dict, use_color_dict = use_color_dict)
    # put in the list
    plot_per_condition[[condition]] <- plot_condition
  }
  return(plot_per_condition)
}


plot_concordance <- function(qtl_output_1, qtl_output_2, snp_column_1, snp_column_2, probe_column_1, probe_column_2, score_column_1, score_column_2, significance_column_1, significance_column_2, assessed_allele_column_1, assessed_allele_column_2){
  # read the tables
  qtls_1 <- fread(qtl_output_1)
  qtls_2 <- fread(qtl_output_2)
  # add columns with standardised names
  qtls_1$SNP <- qtls_1[[snp_column_1]]
  qtls_2$SNP <- qtls_2[[snp_column_2]]
  qtls_1$probe <- qtls_1[[probe_column_1]]
  qtls_2$probe <- qtls_2[[probe_column_2]]
  qtls_1$score_1 <- qtls_1[[score_column_1]]
  qtls_2$score_2 <- qtls_2[[score_column_2]]
  qtls_1$significance_1 <- qtls_1[[significance_column_1]]
  qtls_2$significance_2 <- qtls_2[[significance_column_2]]
  qtls_1$allele_assessed_1 <- qtls_1[[assessed_allele_column_1]]
  qtls_2$allele_assessed_2 <- qtls_2[[assessed_allele_column_2]]
  # subset to what we need, thus speeding it all up
  qtls_1 <- qtls_1[, c('SNP', 'probe', 'allele_assesed_1', 'score_1', 'significance_1'), with=F]
  qtls_2 <- qtls_2[, c('SNP', 'probe', 'allele_assesed_2', 'score_2', 'significance_2'), with=F]
  # set keys for easy joining
  setkey(qtls_1, SNP, probe)
  setkey(qtls_2, SNP, probe)
  qtls <- qtls_1[qtls_2]
  # correct for direction
  qtls[qtls$allele_assesed_1 != qtls$allele_assesed_2, 'score_1', with=F] <- qtls[qtls$allele_assesed_1 != qtls$allele_assesed_2, 'score_1', with=F] * -1
  # TODO continue and plot here
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
  color_coding[["UT_Baseline"]] <- "khaki2"
  color_coding[["UT_t24h"]] <- "khaki4"
  color_coding[["UT_t8w"]] <- "paleturquoise1"
  color_coding[["UT + Baseline"]] <- "khaki2"
  color_coding[["UT + t24h"]] <- "khaki4"
  color_coding[["UT + t8w"]] <- "paleturquoise1"
  color_coding[["Baselinet24h"]] <- "paleturquoise3"
  color_coding[["Baselinet8w"]] <- "rosybrown1"
  color_coding[["t24ht8w"]] <- "rosybrown3"
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["t0-HC"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["t24h-HC"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t8w-HC"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t24h-t0"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t8w-t0"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
  color_coding[["t8w-t24h"]] <- "#C00040" #"rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[['UTBaseline']] <- 'UT-Baseline'
  label_dict[['UTt24h']] <- 'UT-t24h'
  label_dict[['UTt8w']] <- 'UT-t8w'
  label_dict[['Baselinet24h']] <- 'Baseline-t24h'
  label_dict[['Baselinet8w']] <- 'Baseline-t8w'
  label_dict[['t24ht8w']] <- 't24h-t8w'
  label_dict[['UTBaseline']] <- 'HC-t0'
  label_dict[['UTt24h']] <- 'HC-t24h'
  label_dict[['UTt8w']] <- 'HC-t8w'
  label_dict[['Baselinet24h']] <- 't0-t24h'
  label_dict[['Baselinet8w']] <- 't0-t8w'
  # merged combinations
  label_dict[['UT_Baseline']] <- 'UT + Baseline'
  label_dict[['UT_t24h']] <- 'UT + t24h'
  label_dict[['UT_t8w']] <- 'UT + t8w'
  # conditions
  label_dict[['UT']] <- 'HC'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}




# plot locations
egene_sharing_ct_loc <- '/data/cardiology/eQTL_mapping/plots/MatrixEQTL/stemi_and_1mut_lowerres_20210629_eqtlgenlead/egene_sharing_ct/'
egene_sharing_cond_loc <- '/data/cardiology/eQTL_mapping/plots/MatrixEQTL/stemi_and_1mut_lowerres_20210629_eqtlgenlead/egene_sharing_cond/'
# location of output
eqtl_output <- '/data/cardiology/eQTL_mapping/results/MatrixEQTL/stemi_and_1mut_lowerres_20210629_eqtlgenlead/'

# get the condition sharing in each cell type
condition_sharing_per_ct <- plot_egene_sharing_per_condition_allct(eqtl_output, '.cis.fdr.tsv', conditions = c('UT_Baseline', 'UT_t24h', 'UT_t8w'))
# get the cell type sharing in each condition
condition_sharing_per_cond <- plot_egene_sharing_per_celltype_allcond(eqtl_output, conditions = c('UT_Baseline', 'UT_t24h', 'UT_t8w'))

# load object
#cardio.stemi <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.stemi.20210611.combatcorrected.rds')
# genotype location
snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
# load genotype data
snps <- fread(snps_loc, header=T)
# locations of features
features_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/'
# locations of the results
results_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/' 
results_emp_meta_eqtlgenlead <- paste(results_loc, 'inhouse_eQTL_mapping_pipeline/stemi_meta_lowerres_20210301_confine_lead_snp_gene/', sep = '')
#
mash_emp_meta_eqtlgenlead <- get_conditions_cell_types_mash(results_emp_meta_eqtlgenlead)
saveRDS(mash_emp_meta_eqtlgenlead, '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/mashr_stemi_meta_lowerres_20210301_confine_lead_snp_gene.rds')
