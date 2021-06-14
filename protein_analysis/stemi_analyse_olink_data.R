
plot_qtl <- function(expression, genotypes, gene_name, snp_name, paper_style=F){
  # get overlapping expression and genotypes
  participants_both <- intersect(colnames(expression), colnames(expression))
  # get the expression data
  expression_gene <- as.vector(unlist(expression[gene_name, participants_both]))
  # get the snps
  snps <- as.vector(unlist(genotypes[snp_name, participants_both]))
  # turn into plottable dataframe
  plot_table <- data.frame(genotype=snps, expression=expression_gene)
  # turn into a plot
  p <- ggplot(data=plot_table, mapping=aes(x=genotype, y=expression, fill=genotype)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(size = 0.5, alpha = 0.5) +
    ggtitle(paste(snp_name, 'affecting', gene_name))
  # use paper style if requested
  if(paper_style){
    p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  }
  return(p)
}

olink_to_genesymbols <- function(olink, olinkid_to_uid_loc, uniprotid_to_gs_loc){
  # read mapping of olink ID to uniprot ID
  olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T, stringsAsFactors = F)
  # change column names to be uniprot IDs
  colnames(olink) <- olinkid_to_uid[match(colnames(olink), as.character(olinkid_to_uid$OlinkID)), 'Uniprot.ID']
  # read mapping of uniprot IDs to gene symbols
  uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T, stringsAsFactors = F)
  # now change colnames from uniprot IDs to gene symbols
  colnames(olink) <- uniprotid_to_gs[match(colnames(olink), as.character(uniprotid_to_gs$From)), 'To']
  return(olink)
}

olink_to_plottable_table <- function(olink, split_char='\\.'){
  id_and_timepoint_df <- NULL
  # check the rownames
  for(id_and_timepoint in rownames(olink)){
    # split into id and timepoint
    id_and_timepoint_split <- strsplit(id_and_timepoint, split_char)
    id <- id_and_timepoint_split[[1]][[1]]
    timepoint <- id_and_timepoint_split[[1]][[2]]
    # turn into row
    id_and_timepoint_row <- data.frame(id=c(id), timepoint=c(timepoint))
    # add to dataframe
    if(is.null(id_and_timepoint_df)){
      id_and_timepoint_df <- id_and_timepoint_row
    }
    else{
      id_and_timepoint_df <- rbind(id_and_timepoint_df, id_and_timepoint_row)
    }
  }
  # add new columns to original data
  olink <- cbind(id_and_timepoint_df, olink)
  return(olink)
}

plot_olink_expression <- function(olink, identifier, olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL){
  # add extra column
  olink$protein_expression <- olink[[identifier]]
  # set colors based on condition
  colScale <- scale_fill_manual(name = "condition",values = unlist(cc[olink$timepoint]))
  # create the plot
  p <- ggplot(data=olink, mapping=aes(x=timepoint, y=protein_expression, fill=timepoint)) + 
    geom_boxplot(outlier.shape = NA) + 
    colScale +
    geom_jitter(size = 0.5, alpha = 0.5)
  
  title <- paste('protein expression of', identifier)
  if(!is.null(olinkid_to_uid)){
    # read mapping of olink ID to uniprot ID
    olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T, stringsAsFactors = F)
    # get the uid
    uid <- olinkid_to_uid[olinkid_to_uid$OlinkID == identifier, 'Uniprot.ID']
    # add to title
    title <- paste(title, '-', uid, sep = '')
    # if there is a gene symbol mapping, do that one as well
    if(!is.null(uniprotid_to_gs_loc)){
      uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T, stringsAsFactors = F)
      gs <- uniprotid_to_gs[uniprotid_to_gs$From == uid, 'To']
      # add to title
      title <- paste(title, '(', gs, ')')
    }
  }
  # add title
  p <- p + ggtitle(title)
  return(p)
}

plot_olink_expression_per_protein <- function(olink, output_loc, olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL, use_label_dict=F){
  # make labels prettier for timepoint if requested
  if(use_label_dict){
    olink$timepoint <- as.vector(unlist(label_dict()[as.character(olink$timepoint)]))
  }
  # grab the IDs
  identifiers <- setdiff(colnames(olink), c('timepoint', 'id'))
  # plot each protein
  for(identifier in identifiers){
    p <- plot_olink_expression(olink, identifier, olinkid_to_uid_loc, uniprotid_to_gs_loc)
    # set the location
    output_file <- paste(output_loc, identifier, '.pdf', sep = '')
    p
    ggsave(output_file)
  }
}

do_differential_expression_analysis_all <- function(olink, split_char='\\.', method='wilcoxon', paired=T, conditions.x=c('Baseline', 'Baseline'), conditions.y=c('t24h', 't8w'), mtc_method='bonferroni', paired_mtc=T, olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL, inclusion_list_loc=NULL){
  # get it so that the identifier and and timepoint are in the table
  olink <- olink_to_plottable_table(olink, split_char = split_char)
  # include only specific participants if requested
  if(!is.null(inclusion_list_loc)){
    included <- read.table(inclusion_list_loc, header=F, )$V1
    olink <- olink[olink$id %in% included, ]
  }
  # init table
  full_result_table <- NULL
  # check each timepoint
  for(i in 1:length(conditions.x)){
    # grab the conditions
    condition.x <- conditions.x[i]
    condition.y <- conditions.y[i]
    # do the analysis
    result_table <- do_differential_expression_analysis(olink, condition.1 = condition.x, condition.2 = condition.y, timepoint_column='timepoint', sample_column='id',  method=method, paired=paired, olinkid_to_uid_loc=olinkid_to_uid_loc, uniprotid_to_gs_loc=uniprotid_to_gs_loc)
    # to multiple testing if requested
    if(!is.null(mtc_method) & paired_mtc){
      result_table[[paste('p', mtc_method, sep = '_')]] <- p.adjust(result_table$p_nominal, method = mtc_method)
      result_table$paired_mtc <- paired_mtc
    }
    # add to result table
    if(is.null(full_result_table)){
      full_result_table <- result_table
    }
    else{
      full_result_table <- rbind(full_result_table, result_table)
    }
  }
  if(!is.null(mtc_method) & !paired_mtc){
    full_result_table[[paste('p', mtc_method, sep='_')]] <- p.adjust(full_result_table$p_nominal, method = mtc_method)
    full_result_table$paired_mtc <- paired_mtc
  }
  return(full_result_table)
}

do_differential_expression_analysis <- function(olink, condition.1, condition.2, timepoint_column='timepoint', sample_column='id',  method='wilcoxon', paired=T, olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL){
  # grab the samples
  samples_x <- olink[olink[[timepoint_column]] == condition.1, sample_column]
  samples_y <- olink[olink[[timepoint_column]] == condition.2, sample_column]
  # if we go paired, we need to subset to what is in both
  if(paired){
    samples_x <- intersect(samples_x, samples_y)
    samples_y <- intersect(samples_x, samples_y)
  }
  # subset to this data with these samples and that condition
  expression_data_x <- olink[olink[[timepoint_column]] == condition.1 & olink[[sample_column]] %in% samples_x, ]
  expression_data_y <- olink[olink[[timepoint_column]] == condition.2 & olink[[sample_column]] %in% samples_y, ]
  # get the ids of the expression data
  expressed_data <- setdiff(colnames(olink), c(timepoint_column, sample_column))
  # initialize a matrix
  result_matrix <- matrix(, nrow = length(expressed_data), ncol = 10,dimnames = list(expressed_data, c('condition.1', 'condition.2', 'nsamples.x', 'nsamples.y', 'mean_exp.1', 'mean_exp.2', 'mean_lfc', 'method', 'paired', 'p_nominal')))
  # check each expressed datapoint
  for(expressed_datapoint in expressed_data){
    # grab the expression data for x and y
    expression_x <- expression_data_x[match(samples_x, expression_data_x[[sample_column]]), expressed_datapoint]
    expression_y <- expression_data_y[match(samples_y, expression_data_y[[sample_column]]), expressed_datapoint]
    # calculate mean expression
    mean_exp_x <- mean(expression_x)
    mean_exp_y <- mean(expression_y)
    # calculate the mean LFC
    mean_lfc <- log2(mean_exp_x/mean_exp_y)
    # do a statistical analysis
    stat <- NA
    try({
      if(method == 'wilcoxon'){
        res <- wilcox.test(expression_x, expression_y, paired = paired)
        stat <- res$p.value
      }
      else{
        print('unsupported method')
      }
    })
    # add these results
    result_matrix[expressed_datapoint, 'mean_exp.1'] <- mean_exp_x                                      
    result_matrix[expressed_datapoint, 'mean_exp.2'] <- mean_exp_y
    result_matrix[expressed_datapoint, 'mean_lfc'] <- mean_lfc
    result_matrix[expressed_datapoint, 'p_nominal'] <- stat
  }
  # change to dataframe for easier access
  result_matrix <- data.frame(result_matrix)
  if(!is.null(olinkid_to_uid)){
    # read mapping of olink ID to uniprot ID
    olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T, stringsAsFactors = F)
    # get the uid
    uid <- olinkid_to_uid[match(rownames(result_matrix), olinkid_to_uid$OlinkID), 'Uniprot.ID']
    # add to table
    result_matrix$uniprot.id <- uid
    # if there is a gene symbol mapping, do that one as well
    if(!is.null(uniprotid_to_gs)){
      uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T, stringsAsFactors = F)
      gs <- uniprotid_to_gs[match(result_matrix$uniprot.id, uniprotid_to_gs$From), 'To']
      # add to title
      result_matrix$gs <- gs
    }
  }
  # init some data that is the same across the data
  result_matrix$condition.1 <- condition.1
  result_matrix$condition.2 <- condition.2
  result_matrix$nsamples.x <- nrow(expression_data_x)
  result_matrix$nsamples.y <- nrow(expression_data_y)
  result_matrix$method <- method
  result_matrix$paired <- paired
  return(result_matrix)
}

olink_to_emp_format <- function(olink, split_char='\\.', olinkid_to_uid_loc=NULL, uniprotid_to_gs_loc=NULL, gs_to_ens_loc=NULL, inclusion_list_loc=NULL){
  # get it so that the identifier and and timepoint are in the table
  olink <- olink_to_plottable_table(olink, split_char = split_char)
  # include only specific participants if requested
  if(!is.null(inclusion_list_loc)){
    included <- read.table(inclusion_list_loc, header=F, )$V1
    olink <- olink[olink$id %in% included, ]
  }
  # create list with a matrix per condition
  expression_per_condition <- list()
  # check each condition
  for(condition in unique(olink$timepoint)){
    # grab data for that condition
    expression_condition <- olink[olink$timepoint == condition, ]
    # set the sample as row name
    rownames(expression_condition) <- expression_condition$id
    # remove the columns that are not expression data
    expression_condition$id <- NULL
    expression_condition$timepoint <- NULL
    # convert expression names if required
    if(!is.null(olinkid_to_uid_loc)){
      # read mapping of olink ID to uniprot ID
      olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T, stringsAsFactors = F)
      # get the uid
      uid <- olinkid_to_uid[match(colnames(expression_condition), olinkid_to_uid$OlinkID), 'Uniprot.ID']
      # add to table
      colnames(expression_condition) <- uid
      # if there is a gene symbol mapping, do that one as well
      if(!is.null(uniprotid_to_gs_loc)){
        uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T, stringsAsFactors = F)
        gs <- uniprotid_to_gs[match(colnames(expression_condition), uniprotid_to_gs$From), 'To']
        # add to title
        colnames(expression_condition) <- gs
        # if there is ensemble id mapping, to that as well
        if(!is.null(gs_to_ens_loc)){
          gs_to_ens <- read.table(gs_to_ens_loc, sep = '\t')
          ens <- gs_to_ens[match(colnames(expression_condition), gs_to_ens$V2), 'V1']
          colnames(expression_condition) <- ens
        }
      }
    }
    # transpose the table
    expression_condition <- data.frame(t(expression_condition))
    # remove NA
    expression_condition <- expression_condition[!(startsWith(rownames(expression_condition), 'NA.')), ]
    # add to list
    expression_per_condition[[condition]] <- expression_condition
  }
  return(expression_per_condition)
}


correlate_gene_and_protein_expression <- function(protein_expression_matrix, gene_expression_matrix, method='spearman'){
  # get the common participants
  common_parts <- intersect(colnames(protein_expression_matrix), colnames(gene_expression_matrix))
  # get the common genes
  common_genes <- intersect(rownames(protein_expression_matrix), rownames(gene_expression_matrix))
  # check each gene
  correlations <- list()
  for(gene in common_genes){
    # get expression data
    gene_expression <- as.vector(unlist(gene_expression_matrix[gene, common_parts]))
    protein_expression <- as.vector(unlist(protein_expression_matrix[gene, common_parts]))
    # calculate the correlation
    correlation <- NA
    try({
      # inside a try block so if we can't correlate, we'll still do the rest
      correlation <- cor(gene_expression, protein_expression)
    })
    correlations[[gene]] <- correlation
  }
  return(correlations)
}


correlate_gene_and_protein_expression_per_cell_type <- function(protein_expression_location, gene_expression_location, cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  # create a dataframe to store results
  correlation_per_cell_type <- NULL
  # load the protein expression
  protein_expression <- read.table(protein_expression_location, sep = '\t', header = T, row.names = 1)
  # check each cell type
  for(cell_type in cell_types){
    # make the location of the cell type
    gene_expression_cell_type_loc <- paste(gene_expression_location, cell_type, '_expression.tsv', sep = '')
    # read that table
    gene_expression <- read.table(gene_expression_cell_type_loc, header = T, row.names = 1)
    # get the correlations
    correlations <- correlate_gene_and_protein_expression(protein_expression, gene_expression)
    # turn into a dataframe
    correlation_columns <- data.frame(gene=names(correlations), correlation=as.vector(unlist(correlations)))
    # set the correlations as cell type
    colnames(correlation_columns) <- c('gene', cell_type)
    # add to existing dataframe
    if(is.null(correlation_per_cell_type)){
      correlation_per_cell_type <- correlation_columns
    }
    else{
      correlation_per_cell_type <- merge(correlation_per_cell_type, correlation_columns, by='gene', all=T)
    }
  }
  return(correlation_per_cell_type)
}


correlate_gene_and_protein_expression_per_condition <- function(protein_expression_location, gene_expression_location, conditions=c('Baseline', 't24h', 't8w'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  # save table per condition
  table_per_condition <- list()
  # check each condition
  for(condition in conditions){
    # paste the protein location together
    protein_expression_location_condition <- paste(protein_expression_location, condition, '/bulk_expression.tsv', sep = '')
    # paste the gene expression location together
    gene_expression_location_condition <- paste(gene_expression_location, condition, '/', sep = '')
    # save result in list
    table_per_condition[[condition]] <- correlate_gene_and_protein_expression_per_cell_type(protein_expression_location, gene_expression_location_condition, cell_types)
  }
  return(table_per_condition)
}

get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UTBaseline"]] <- "khaki2"
  color_coding[["UTt24h"]] <- "khaki4"
  color_coding[["UTt8w"]] <- "paleturquoise1"
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
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "paleturquoise3"
  color_coding[["t0-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["Baseline"]] <- "pink"
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


olinkid_to_uid_loc <- '/data/cardiology/olinkid_to_uniprotid.tsv'
uniprotid_to_gs_loc <- '/data/cardiology/uniprot_to_genesymbol.tsv'
gs_to_ens_loc <- '/data/cardiology/eQTL_mapping/features_v3.tsv'
olink_loc <- '/data/cardiology/20200442_Groot_NPX-QC_format_fixed.tsv'
features_loc <- '/data/cardiology/olink/features/stemi_features_20210608/'
inclusion_list_loc <- '/data/cardiology/included_participants.txt'

# read with rownames instead
olink <- read.table(olink_loc, sep = '\t', header = T, row.names = 1)
olink_plottable <- olink_to_plottable_table(olink)
# remove the test samples
olink_plottable <- olink_plottable[olink_plottable$timepoint %in% c('Baseline', 't24h', 't8w'), ]
plot_olink_expression_per_protein(olink_plottable, '/data/cardiology/olink/plots/expression/', olinkid_to_uid_loc, uniprotid_to_gs_loc)
# do DE analysis
olink_de <- do_differential_expression_analysis_all(olink=olink, olinkid_to_uid_loc = olinkid_to_uid_loc, uniprotid_to_gs_loc = uniprotid_to_gs_loc, inclusion_list_loc=inclusion_list_loc)
write.table(olink_de, '/data/cardiology/differential_expression/olink_de_paired_wilcoxon_baseline.tsv', row.names=T, col.names = T, sep = '\t')
# convert to EMP style tables
olink_emp <- olink_to_emp_format(olink, olinkid_to_uid_loc = olinkid_to_uid_loc, uniprotid_to_gs_loc = uniprotid_to_gs_loc, gs_to_ens_loc = gs_to_ens_loc, inclusion_list_loc=inclusion_list_loc)
write.table(olink_emp[['Baseline']], paste(features_loc, 'Baseline/bulk_expression.tsv', sep=''), sep = '\t', row.names=T, col.names = NA, quote = F)
write.table(olink_emp[['t24h']], paste(features_loc, 't24h/bulk_expression.tsv', sep=''), sep = '\t', row.names=T, col.names = NA, quote = F)
write.table(olink_emp[['t8w']], paste(features_loc, 't8w/bulk_expression.tsv', sep=''), sep = '\t', row.names=T, col.names = NA, quote = F)
