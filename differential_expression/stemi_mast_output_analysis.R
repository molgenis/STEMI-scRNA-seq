####################
# libraries        #
####################

library(RColorBrewer)
library(MetaVolcanoR)
library(ggplot2)
library(data.table)
library(ggpubr)
library(UpSetR)
library(VennDiagram)
library(familyR)
library(visNetwork)

####################
# Functions        #
####################

write_version_chemistry_overlap <- function(mast_meta_output_loc, venn_output_loc){
  # grab colors
  myCol <- brewer.pal(3, "Pastel2")
  # check each meta DE MAST file
  files <- list.files(mast_meta_output_loc)
  for(file in files){
    # read the file
    table <- read.table(paste(mast_meta_output_loc, '/', file, sep = ''), sep = '\t', header=T)
    # create an output name for the venn diagram we'll create
    output_file <- substr(file, 1, regexpr("\\.[^\\.]*$", file)-1)
    # check for overlap of significant DE genes in v2, v3 and meta
    venn.diagram(x = list(rownames(table[table$p_val_v2*nrow(table)<0.05, ]),rownames(table[table$p_val_v3*nrow(table)<0.05, ]),rownames(table[table$metap_bonferroni<0.05, ])),
                 category.names = c('v2', 'v3', 'meta'),
                 main = substr(file, 1, regexpr("\\.[^\\.]*$", file)-1),
                 filename = paste(venn_output_loc, '/', output_file, '.png', sep = ''),
                 imagetype="png" ,
                 height = 480 , 
                 width = 480 , 
                 resolution = 300,
                 compression = "lzw",
                 lwd = 2,
                 lty = 'blank',
                 fill = myCol,
                 cex = .6,
                 fontface = "bold",
                 fontfamily = "sans",
                 cat.cex = 0.6,
                 cat.fontface = "bold",
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27, 135),
                 cat.dist = c(0.055, 0.055, 0.085),
                 cat.fontfamily = "sans",
                 rotation = 1)
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
        # get the meta p values using stouffers method
        #stouffers <- rep(NA, times = nrow(mast))
        #for(i in 1:nrow(mast)){
        #  # get the p-values
        #  p_vals <- c(mast[i, 'p_val_v2'], mast[i, 'p_val_v3'])
        #  # the weights are based on the number of cells
        #  weights <- c(sqrt(cond1_v2_cells + cond2_v2_cells), sqrt(cond1_v3_cells + cond2_v3_cells))
        #  # get the result from the Stouffer's method
        #  stouffers_res <- sumz(p = p_vals, weights = weights)
        #  if(!is.na(stouffers_res)){
        #    stouffers[i] <- stouffers_res$p[1,1]*length(genes_both) #bonferroni correct by multiplying by number of tests
        #  }
        #}
        # add the value
        #mast$stouffers_p <- stouffers
        #mast[mast$stouffers_p > 1 & !is.na(mast$stouffers_p), ]$stouffers_p <- 1
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


get_significant_genes <- function(mast_output_loc, sig_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # get the files
  files <- list.files(mast_output_loc)
  # try to read each file
  for(file in files){
    try({
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # confine in some way if reporting a max number of genes
      if(!is.null(max)){
        # by p if required
        if(max_by_pval){
          mast <- mast[order(mast[[p_val_column]]), ]
        }
        # by lfc otherwise
        else{
          mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        }
        # subset to the number we requested if max was set
        mast <- mast[1:max,]
      }
      # grab the genes from the column names
      genes <- rownames(mast)
      # convert the symbols to ensemble IDs
      if (to_ens) {
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        genes <- mapping[match(genes, mapping$V2),"V1"]
      }
      # otherwise change the Seurat replacement back
      else{
        #genes <- gsub("-", "_", genes)
      }
      # create a regex to get the last index of .
      last_dot_pos <- "\\.[^\\.]*$"
      # this allows us to remove the filename extention
      file_no_ext <- substr(file, 1, regexpr(last_dot_pos,file)-1)
      # create output location
      sig_output <- paste(sig_output_loc, file_no_ext, '.txt', sep = '')
      # write the genes
      write.table(genes, sig_output, sep = '\t', quote = F, row.names = F, col.names = F)
    })
  }
}

get_pathway_table <- function(pathway_output_loc, sig_val_to_use = 'q.value.FDR.B.H', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), append='_sig_pathways.txt', use_ranking=T){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- paste(pathway_output_loc, cell_type, stim,stim2, append, sep = '')
            # read the file
            pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
            # create column name
            newcolname <- paste(cell_type, stim, stim2, sep = '')
            # use ranking or log value
            if(use_ranking){
              pathways[[newcolname]] <- as.numeric(rownames(pathways))
            }
            else{
              pathways[[newcolname]] <- (pathways[[sig_val_to_use]])
            }
            #pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
            pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
            # reduce to the only two columns we care about
            pathways <- pathways[, c('id_name', newcolname)]
            # join with other pathway files
            if(is.null(pathway_df)){
              # just set as df if the first round through
              pathway_df <- pathways
              pathway_df <- data.table(pathway_df, key = c('id_name'))
            }
            else{
              # otherwise, merge with existing pathways
              pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
              #pathway_df[[newcolname]] <- pathways[[newcolname]][match(pathway_df$Name, pathways$Name)]
              #pathway_df <- left_join(pathway_df, pathways)
              
            }
          }
        })
      }
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}

get_pathways <- function(pathway_output_loc, append='_sig_pathways.txt', sig_val_to_use = 'q.value.FDR.B.H', sig_pval=0.05, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w')){
  # set up per cell type
  pathways_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # set up per stim combination
    pathways_per_condition_combination <- list()
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- paste(pathway_output_loc, cell_type, stim,stim2, append, sep = '')
            # read the file
            pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
            # filter to only the significant rows
            pathways <- pathways[pathways[[sig_val_to_use]] < sig_pval, ]
            # grab the combined names
            pathway_names <- paste(pathways$ID, pathways$Name, sep = '_')
            # put in the list
            pathways_per_condition_combination[[paste(stim, stim2, sep='')]] <- pathway_names
          }
        })
      }
    }
    # put in the list per cell type
    pathways_per_ct[[cell_type]] <- pathways_per_condition_combination
  }
  return(pathways_per_ct)
}


pathway_numbers_to_table <- function(pathway_output_loc, append='_sig_pathways.txt', sig_val_to_use = 'q.value.FDR.B.H', sig_pval=0.05, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), remove_na_cols=T){
  # get the DE genes
  pathways_per_ct <- get_pathways(pathway_output_loc=pathway_output_loc, append=append, sig_val_to_use=sig_val_to_use, sig_pval=sig_pval, cell_types=cell_types, stims=stims)
  # make all possible combinations of stims
  combs <- paste(rep(stims, each = length(stims)), stims, sep = '')
  # create matrix to store results
  number_table <- matrix(, ncol=length(combs), nrow=length(cell_types), dimnames = list(cell_types, combs))
  # check each cell type
  for(cell_type in intersect(cell_types, names(pathways_per_ct))){
    # get for specific cell type
    pathways_per_conditin <- pathways_per_ct[[cell_type]]
    # check each condition combination
    for(comb in intersect(combs, names(pathways_per_conditin))){
      # get the number of genes
      nr_of_pathways <- length(pathways_per_conditin[[comb]])
      # add to matrix
      number_table[cell_type, comb] <- nr_of_pathways
    }
  }
  # remove na column
  if(remove_na_cols){
    number_table <- number_table[, colSums(is.na(number_table)) < nrow(number_table)]
  }
  # i like dataframes
  number_table <- data.frame(number_table)
  return(number_table)
}

get_de_genes <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w')){
  # set up per cell type
  de_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # set up per stim combination
    de_per_condition <- list()
    # check each stim
    for(stim in stims){
      for(stim2 in stims){
        try({
          if(stim != stim2){
            print(paste(cell_type, stim, stim2, sep = ' '))
            # paste the filepath together
            filepath <- paste(mast_output_loc, cell_type, stim,stim2, '.tsv', sep = '')
            # read the file
            # read the mast output
            mast <- read.table(filepath, header=T)
            # filter to only include the significant results
            mast <- mast[mast[[pval_column]] <= 0.05, ]
            # filter for only the positive lfc if required
            if(only_positive){
              mast <- mast[mast[[lfc_column]] < 0, ]
            }
            # filter for only the positive lfc if required
            if(only_negative){
              mast <- mast[mast[[lfc_column]] > 0, ]
            }
            # confine in some way if reporting a max number of genes
            if(!is.null(max)){
              # by p if required
              if(max_by_pval){
                mast <- mast[order(mast[[p_val_column]]), ]
              }
              # by lfc otherwise
              else{
                mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
              }
              # subset to the number we requested if max was set
              mast <- mast[1:max,]
            }
            # grab the genes from the column names
            genes <- rownames(mast)
            # convert the symbols to ensemble IDs
            if (to_ens) {
              mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
              mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
              genes <- mapping[match(genes, mapping$V2),"V1"]
            }
            # otherwise change the Seurat replacement back
            else{
              #genes <- gsub("-", "_", genes)
            }
            de_per_condition[[paste(stim, stim2, sep = '')]] <- genes
          }
        })
      }
    }
    de_per_ct[[cell_type]] <- de_per_condition
  }
  return(de_per_ct)
}


de_genes_number_to_table <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w'), remove_na_cols=T){
  # get the DE genes
  de_genes_per_ct <- get_de_genes(mast_output_loc, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping, cell_types=cell_types, stims=stims)
  # make all possible combinations of stims
  combs <- paste(rep(stims, each = length(stims)), stims, sep = '')
  # create matrix to store results
  number_table <- matrix(, ncol=length(combs), nrow=length(cell_types), dimnames = list(cell_types, combs))
  # check each cell type
  for(cell_type in intersect(cell_types, names(de_genes_per_ct))){
    # get for specific cell type
    de_genes_per_conditin <- de_genes_per_ct[[cell_type]]
    # check each condition combination
    for(comb in intersect(combs, names(de_genes_per_conditin))){
      # get the number of genes
      nr_of_de_genes <- length(de_genes_per_conditin[[comb]])
      # add to matrix
      number_table[cell_type, comb] <- nr_of_de_genes
    }
  }
  # remove na column
  if(remove_na_cols){
    number_table <- number_table[, colSums(is.na(number_table)) < nrow(number_table)]
  }
  # i like dataframes
  number_table <- data.frame(number_table)
  return(number_table)
}

numbers_table_to_plot <- function(numbers_table, cols_include=NULL, use_label_dict=T, use_groups_dict=T, title=NULL, pointless=F, legendless=F){
  numbers_table_to_do <- numbers_table
  if(!is.null(cols_include)){
    numbers_table_to_do <- numbers_table_to_do[, cols_include, drop = F]
  }
  plot_data <- NULL
  for(cell_type in rownames(numbers_table_to_do)){
    for(condition_comb in colnames(numbers_table_to_do)){
      val <- numbers_table_to_do[cell_type, condition_comb]
      if(is.na(val)){
        val <- 0
      }
      # set labels to use in plot
      conditions_label <- condition_comb
      cell_type_label <- cell_type
      # create nicer labels if requested
      if(use_label_dict){
        conditions_label <- label_dict()[[conditions_label]]
        cell_type_label <- label_dict()[[cell_type]]
      }
      # create the row
      row_plot_data <- data.frame(de_genes=val, cell_type=cell_type_label, conditions=conditions_label)
      if(use_groups_dict){
        row_plot_data$comparison <- groups_dict()[[condition_comb]]
      }
      if(is.null(plot_data)){
        plot_data <- row_plot_data
      }
      else{
        plot_data <- rbind(plot_data, row_plot_data)
      }
    }
  }
  # create the ylim
  ylims <- c(0, max(plot_data$de_genes*1.1))
  # make the plots
  p <- NULL
  if(use_groups_dict){
    p <- ggplot(data=plot_data, aes(x=comparison, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type) + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) + ylim(ylims)
  }
  else if(length(unique(plot_data$conditions))==1){
    p <- ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=cell_type), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(values = unlist(get_color_coding_dict()[unique(plot_data$cell_type)])) +
      xlab('') + 
      ylab('number of DE genes') +
      theme(axis.text.x=element_blank(), 
            axis.ticks = element_blank(), 
            legend.title = element_text(size=14), 
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=12),
            strip.text.x = element_text(size=12),
            legend.position = 'none') + ylim(ylims)
  }
  else{
    #ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) + theme(legend.position = 'none') +
    p <- ggplot(data=plot_data, aes(x=conditions, y=de_genes)) + geom_point(aes(color=conditions), size=6) + facet_grid(. ~ cell_type)  + scale_color_manual(name = 'condition\ncombination', values = unlist(get_color_coding_dict()[unique(plot_data$conditions)])) +
      xlab('condition combination') + 
      ylab('number of DE genes') +
      theme(#axis.text.x=element_blank(), 
            #axis.ticks = element_blank(), 
            legend.title = element_text(size=14), 
            legend.text = element_text(size=12),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14),
            axis.text.y = element_text(size=12),
            strip.text.x = element_text(size=12)) + ylim(ylims)
  }
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }
  if(pointless){
    p <- p + theme(axis.text.x=element_blank(), 
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank())
  }
  if(legendless){
    p <- p + theme(legend.position = 'none')
  }
  return(p)
}

split_label_dict <- function(){
  label_dict <- list()
  #label_dict[['UTBaseline']] <- 'UT\nBaseline'
  #label_dict[['UTt24h']] <- 'UT\nt24h'
  #label_dict[['UTt8w']] <- 'UT\nt8w'
  #label_dict[['Baselinet24h']] <- 'Baseline\nt24h'
  #label_dict[['Baselinet8w']] <- 'Baseline\nt8w'
  #label_dict[['t24ht8w']] <- 't24h\nt8w'
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
  return(label_dict)
}

groups_dict <- function(){
  groups_dict <- list()
  groups_dict[['UTBaseline']] <- 'HC'
  groups_dict[['UTt24h']] <- 'HC'
  groups_dict[['UTt8w']] <- 'HC'
  groups_dict[['Baselinet24h']] <- 'STEMI'
  groups_dict[['Baselinet8w']] <- 'STEMI'
  groups_dict[['t24ht8w']] <- 'STEMI'
  return(groups_dict)
}

get_top_pathways <- function(pathway_table, nr_of_top_genes, is_ranked=F){
  # init pathways list
  pathways <- c()
  # go through the columns
  for(col in colnames(pathway_table)){
    # order by that column
    ordered <- pathway_table[order(pathway_table[[col]], decreasing = T), ]
    if(is_ranked){
      ordered <- pathway_table[order(pathway_table[[col]], decreasing = F), ]
    }
    # get those top ones
    top_col <- rownames(ordered)[1:nr_of_top_genes]
    pathways <- c(pathways, top_col)
  }
  # limit to those top pathways now
  pathway_table_smaller <- pathway_table[rownames(pathway_table) %in% pathways, ]
  return(pathway_table_smaller)
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
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
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


get_top_de_genes_per_cond_and_ct <- function(mast_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=5, max_by_pval=T, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', 'Baseline', 't24h', 't8w')){
  # create dataframe for cell type
  df_top_genes <- NULL
  # make all possible combinations of stims
  combs <- paste(rep(stims, each = length(stims)), stims, sep = '')
  # check each cell type
  for(cell_type in cell_types){
    # check each combination
    for(comb in combs){
      # paste the output location together
      output_loc_full <- paste(mast_output_loc, cell_type, comb, '.tsv', sep = '')
      # get the up genes
      up_de_genes <- get_de_genes_from_mast_file(mast_full_file_path = output_loc_full, pval_column = pval_column, sig_pval = sig_pval, max = max, max_by_pval = max_by_pval, lfc_column = lfc_column, to_ens = to_ens, symbols.to.ensg.mapping = symbols.to.ensg.mapping, only_positive = T)
      down_de_genes <- get_de_genes_from_mast_file(mast_full_file_path = output_loc_full, pval_column = pval_column, sig_pval = sig_pval, max = max, max_by_pval = max_by_pval, lfc_column = lfc_column, to_ens = to_ens, symbols.to.ensg.mapping = symbols.to.ensg.mapping, only_negative = T)
      # add to top genes if possible
      if(!is.null(up_de_genes)){
        # turn into df
        up_de_genes_df <- data.frame(a = up_de_genes)
        colnames(up_de_genes_df) <- paste(cell_type, comb, 'up', sep='_')
        if(is.null(df_top_genes)){
          df_top_genes <- up_de_genes_df
        }
        else{
          df_top_genes <- cbind(df_top_genes, up_de_genes_df)
        }
      }
      if(!is.null(down_de_genes)){
        # turn into df
        down_de_genes_df <- data.frame(a = down_de_genes)
        colnames(down_de_genes_df) <- paste(cell_type, comb, 'down', sep='_')
        if(is.null(df_top_genes)){
          df_top_genes <- down_de_genes_df
        }
        else{
          df_top_genes <- cbind(df_top_genes, down_de_genes_df)
        }
      }
    }
  }
  return(df_top_genes)
}


get_top_vary_genes <- function(de_table, use_tp=T, use_ct=T, sd_cutoff=0.5, use_dynamic_sd=F, top_so_many=10, must_be_positive_once=F, timepoints=c("Baselinet24h", "Baselinet8w", "t24ht8w"), cell_types=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")){
  top_vary_de <- c()
  cols_to_loop <- NULL
  # grab the appriate grep
  if(use_tp & use_ct){
    # cell type at a timepoint, so monocyte3hCA and monocyte3hPA and monocyte3hMTB for example
    cols_to_loop <- paste(rep(cell_types, each = length(timepoints)), timepoints, sep = ".*")
  }
  else if(use_tp){
    cols_to_loop <- timepoints
  }
  else if(use_ct){
    cols_to_loop <- cell_types
  }
  else{
    cols_to_loop <- c(paste(colnames(de_table), collapse='|'))
  }
  # go through our group of columns
  for(col_grep in cols_to_loop){
    # grab the column names that have this in their name
    appropriate_columns <- colnames(de_table)[(grep(col_grep, colnames(de_table)))]
    print('getting most varying out of: ')
    print(appropriate_columns)
    # now subset the frame to only have these columns
    sub_de_table <- de_table[, appropriate_columns]
    # subset to only the genes that were upregulated at least once, if requested
    if(must_be_positive_once){
      sub_de_table <- sub_de_table[apply(sub_de_table,1,min) < 0,]
    }
    # we will return the rownames
    varying_genes <- NULL
    # either use a set SD or grab so many genes
    if(use_dynamic_sd){
      varying_genes <- get_most_varying_from_df(sub_de_table, top_so_many)
    }
    else{
      # now calculate the sd over this set of columns
      sds <- apply(sub_de_table, 1, sd, na.rm=T)
      # then grab the genes that are 'this' varied
      varying_genes <- rownames(sub_de_table[sds > sd_cutoff,])
    }
    
    # and add them to the list
    top_vary_de <- c(top_vary_de, varying_genes)
  }
  # constrain to the unique genes
  top_vary_de <- unique(top_vary_de)
  top_vary_de <- sort(top_vary_de)
  return(top_vary_de)
}

get_most_varying_from_df <- function(dataframe, top_so_many=10){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = T), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
}

get_combined_meta_de_table <- function(meta_output_loc, must_be_positive_once=F, timepoints=c("Baselinet24h", "Baselinet8w", "t24ht8w"), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")){
  de_table <- read.table(paste0(meta_output_loc, cell_types_to_use[1], timepoints[1], ".tsv"), stringsAsFactors = F, sep = "\t")
  de_table <- de_table['metafc']
  colnames(de_table) <- c(paste(cell_types_to_use[1], timepoints[1], sep=''))
  deg_meta_combined <- data.frame(row.names = rownames(de_table))
  rows <- rownames(deg_meta_combined)
  deg_meta_combined <- data.table(deg_meta_combined)
  deg_meta_combined$genes <- rows
  
  
  for (timepoint in timepoints) {
    for (cell_type in cell_types_to_use) {
      deg_table <- read.table(paste0(meta_output_loc, cell_type, timepoint, ".tsv"), stringsAsFactors = F, sep = "\t")
      deg_table <- deg_table['metafc']
      colnames(deg_table) <- c(paste(cell_type, timepoint, sep=''))
      deg_table$genes <- rownames(deg_table)
      deg_table <- data.table(deg_table)
      print(head(deg_table))
      deg_meta_combined <- merge(deg_meta_combined, deg_table, by.x='genes', by.y='genes', all=TRUE)
      #deg_meta_combined[,paste(cell_type, timepoint, pathogen, sep = "_")] <- deg_table$metafc
    }
  }
  
  deg_meta_combined <- data.frame(deg_meta_combined)
  rownames(deg_meta_combined) <- deg_meta_combined$genes
  deg_meta_combined$genes <- NULL
  deg_meta_combined[is.na(deg_meta_combined)] <- 0
  # limit to those that were upregulated at least once if requested
  if(must_be_positive_once){
    deg_meta_combined <- deg_meta_combined[apply(deg_meta_combined,1,min) < 0,]
  }
  return(deg_meta_combined)
}


plot_de_vs_cell_type_numbers <- function(mast_output_loc, metadata, timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', 'Baselinet24h', 'Baselinet8w', 't24ht8w'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, plot_separately=F, proportion=F, use_label_dict=T){
  # init table
  table <- NULL
  for(timepoint in timepoints){
    for(cell_type in cell_types_to_use){
      # grab the number of cells, I would like to do this in a better way, but can't think of something now
      nr_of_cells <- 0
      if(timepoint == 'UTBaseline'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 'UT', 'Baseline', proportion = proportion)
      }
      else if(timepoint == 'UTt24h'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 'UT', 't24h', proportion = proportion)
      }
      else if(timepoint == 'UTt8w'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 'UT', 't8w', proportion = proportion)
      }
      else if(timepoint == 'Baselinet24h'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 'Baseline', 't24h', proportion = proportion)
      }
      else if(timepoint == 'Baselinet8w'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 'Baseline', 't8w', proportion = proportion)
      }
      else if(timepoint == 't24ht8w'){
        nr_of_cells <- get_nr_of_cells(metadata, cell_type, 't24h', 't8w', proportion = proportion)
      }
      # grab the nr of genes that are significant
      de_table_loc <- paste(mast_output_loc, cell_type, timepoint, '.tsv', sep = '')
      de_table <- read.table(de_table_loc, sep = '\t', header = T, row.names = 1)
      sig_gene_number <- nrow(de_table[de_table[[pval_column]] < sig_pval, ])
      # add to the table
      if(is.null(table)){
        table <- data.frame(c(timepoint), c(cell_type), c(nr_of_cells), c(sig_gene_number), stringsAsFactors = F)
        colnames(table) <- c('timepoint', 'cell_type', 'nr_of_cells', 'sig_gene_number')
      }
      else{
        table <- rbind(table, c(timepoint, cell_type, nr_of_cells, sig_gene_number))
      }
    }
    
    if(plot_separately){
      ggplot(table[table$timepoint==timepoint, ], aes(x=nr_of_cells, y=sig_gene_number, shape=timepoint, color=cell_type)) +
        geom_point() +
        labs(x = xlab, y = 'nr of DE genes', title = 'cells vs nr of DE genes', shape='condition combination', color='cell type')
    }
  }
  print(table)
  if(!plot_separately){
    # the label is different depending on whether we show the proportion or not
    xlab <- 'nr of cells'
    if(proportion){
      xlab <- 'proportion of cells'
    }
    if(use_label_dict){
      label_dict <- split_label_dict()
      table$timepoint <- as.vector(unlist(label_dict[table$timepoint]))
    }
    # get color coding
    cc <- get_color_coding_dict()
    colScale <- scale_colour_manual(name = "cell type",values = unlist(cc[cell_types_to_use]))
    ggplot(table, aes(x=as.numeric(nr_of_cells), y=as.numeric(sig_gene_number), shape=timepoint, color=cell_type)) +
      #scale_x_discrete(breaks = seq(0, 1, by = 0.1)) +
      #scale_y_discrete(breaks = seq(0, 1000, by = 100)) +
      geom_point(size=3) +
      labs(x = xlab, y = 'nr of DE genes', title = 'cells vs nr of DE genes', shape='condition combination') +
      colScale
  }
}

get_nr_of_cells <- function(metadata, cell_type, condition1, condition2, cell_type_column='cell_type_lowerres', condition_column='timepoint.final', proportion=F){
  # the number of cells is the number of rows in the metadata that are of that cell type
  nr_of_cells <- nrow(metadata[metadata[[cell_type_column]] == cell_type & (metadata[[condition_column]] == condition1 | metadata[[condition_column]] == condition2), ])
  # if we want the proportion, we need to divide by the total number of cells, which is all rows
  if(proportion){
    nr_of_cells <- nr_of_cells / nrow(metadata[metadata[[condition_column]] == condition1 | metadata[[condition_column]] == condition2, ])
  }
  return(nr_of_cells)
}

filter_pathway_df_on_starting_id <- function(pathway_df, filtered_pathway_names, remove_id_from_pathway_name=T){
  # get the ones now in the pathway df
  pathway_names_with_id <- rownames(pathway_df)
  last_dash_pos <- "\\_[^\\_]*$"
  # get the pathway names by skipping from the underscore
  pathway_names <- substr(pathway_names_with_id, regexpr(last_dash_pos, pathway_names_with_id)+1, nchar(pathway_names_with_id))
  print(head(pathway_names))
  # filter the pathway df
  pathway_df_filtered <- pathway_df[pathway_names %in% filtered_pathway_names, ]
  # remove the ID from the pathway name if asked
  if(remove_id_from_pathway_name){
    rownames(pathway_df_filtered) <- substr(rownames(pathway_df_filtered), regexpr(last_dash_pos, rownames(pathway_df_filtered))+1, nchar(rownames(pathway_df_filtered)))
  }
  return(pathway_df_filtered)
}

get_filtered_pathway_names <- function(pathway_table, relation_table, starting_id){
  # get all of the children of the starting ID
  all_children <- get_children(relation_table, starting_id)
  # get the names of the pathways that are children
  pathway_names <- as.character(pathway_table[pathway_table$V1 %in% all_children, ]$V2)
  return(pathway_names)
}

get_children <- function(relation_table, starting_id){
  # get all of the children of the starting ID
  children <- as.character(relation_table[relation_table$V1 == starting_id, 'V2'])
  # these children are all family
  family <- children
  # see if there were any children
  if(length(children) > 0){
    # if there were children, we need to get their children as well
    for(child in children){
      # get the grandchildren and add these to the family
      grand_children <- get_children(relation_table, child)
      family <- c(family, grand_children)
    }
  }
  return(family)
}

get_subcell_ratio <- function(cell_type_large, subcell_type, metadata, condition, cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', condition.column='timepoint.final'){
  number_of_cells_large <- nrow(metadata[metadata[[condition.column]] == condition & metadata[[cell_type_column_lowerres]] == cell_type_large, ])
  number_of_cells_small <- nrow(metadata[metadata[[condition.column]] == condition & metadata[[cell_type_column_higherres]] == subcell_type, ])
  ratio <- number_of_cells_small/number_of_cells_large
  return(ratio)
}


plot_de_number_vs_subcell_population <- function(mast_output_loc, cell_type_large, subcell_type, metadata, timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', condition.column='timepoint.final', pval_column='metap_bonferroni', sig_pval=0.05, make_absolute=F){
  # init table
  table <- NULL
  for(timepoint in timepoints){
      # grab the number of cells, I would like to do this in a better way, but can't think of something now
      ratio1 <- NA
      ratio2 <- NA
      if(timepoint == 'UTBaseline'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'UT', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'Baseline', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      else if(timepoint == 'UTt24h'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'UT', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't24h', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      else if(timepoint == 'UTt8w'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'UT', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't8w', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      else if(timepoint == 'Baselinet24h'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'Baseline', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't24h', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      else if(timepoint == 'Baselinet8w'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 'Baseline', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't8w', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      else if(timepoint == 't24ht8w'){
        ratio1 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't24h', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
        ratio2 <- get_subcell_ratio(cell_type_large, subcell_type, metadata, 't8w', cell_type_column_higherres, cell_type_column_lowerres, condition.column)
      }
      # switch to log2 so it is a positive/negative number instead of 0-1 and 1+
      log2ratio1 <- log2(ratio1)
      log2ratio2 <- log2(ratio2)
      # calculate difference
      difference <- diff(c(log2ratio1,log2ratio2))
      # convert to absolute number if requested
      if(make_absolute){
        difference <- abs(difference)
      }
      # grab the nr of genes that are significant
      de_table_loc <- paste(mast_output_loc, cell_type_large, timepoint, '.tsv', sep = '')
      de_table <- read.table(de_table_loc, sep = '\t', header = T, row.names = 1)
      sig_gene_number <- nrow(de_table[de_table[[pval_column]] < sig_pval, ])
      # add to the table
      if(is.null(table)){
        table <- data.frame(c(timepoint), c(cell_type_large), c(difference), c(sig_gene_number), stringsAsFactors = F)
        colnames(table) <- c('timepoint', 'cell_type', 'log2_ratio_difference', 'sig_gene_number')
      }
      else{
        table <- rbind(table, c(timepoint, cell_type_large, difference, sig_gene_number))
      }
  }
  ggplot(table, aes(x=as.numeric(log2_ratio_difference), y=as.numeric(sig_gene_number), color=timepoint)) +
    geom_point(size=3) +
    labs(x = 'log2 subceltype ratio difference', y = 'nr of DE genes', title = 'subcelltype ratio difference vs nr of DE genes')   
}


plot_DE_sharing_per_celltype <- function(condition_combination, mast_output_loc, cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, only_positive=F, only_negative=F, lfc_column='metafc', use_label_dict=T, use_color_dict=T){
  DE_genes_per_ct <- list()
  # get the DE genes for each cell type
  for(cell_type in cell_types_to_use){
    # build the full path
    full_mast_path <- paste(mast_output_loc, cell_type, condition_combination, '.tsv', sep = '')
    # grab the significant genes
    try({
      # read the mast output
      mast <- read.table(full_mast_path, header=T, row.names = 1, sep = '\t')
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # we just care about the gene names
      sig_genes <- rownames(mast)
      # store these for the cell type
      DE_genes_per_ct[[cell_type]] <- sig_genes
    })
  }
  if(use_label_dict){
    names(DE_genes_per_ct) <- label_dict()[names(DE_genes_per_ct)]
  }
  queries <- list()
  sets.bar.color <- 'black'
  if(use_color_dict){
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    for(i in 1:length(names(DE_genes_per_ct))){
      cell_type <- names(DE_genes_per_ct)[i]
      # add for the singles in the intersection sizes
      ct_list <- list(
        query = intersects,
        params = list(cell_type),
        color = get_color_coding_dict()[[cell_type]],
        active = T)
      queries[[i]] <- ct_list
      # add for the DF to order the set sizes
      numbers_row <- data.frame(ct=c(cell_type), nr=c(length(DE_genes_per_ct[[cell_type]])), stringsAsFactors = F)
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
  upset(fromList(DE_genes_per_ct), order.by = 'freq', nsets = length(DE_genes_per_ct), queries = queries, sets.bar.color=sets.bar.color	)
  #return(DE_genes_per_ct)
}

pathway_mapping_filtered_childless <- function(list_of_pathways_and_parents, pathway_mapping){
  # full list of terms
  all_terms <- c()
  for(start in names(list_of_pathways_and_parents)){
    # get all the terms
    all_terms <- c(all_terms, start, list_of_pathways_and_parents[[start]])
  }
  # make the unique terms
  all_terms <- unique(all_terms)
  # iteratively remove parents with no children in our data
  keep_going <- T
  if(keep_going == T){
    # get the childless entries
    child_less_parents <- setdiff(pathway_mapping$V2, pathway_mapping$V1)
    # check if any of those are not in our set
    child_less_parents_filtered <- setdiff(child_less_parents, all_terms)
    # if there are none left, stop working
    if(length(child_less_parents_filtered) == 0){
      keep_going <- F
    }
    else{
      # remove these children without children that we don't care about
      pathway_mapping <- pathway_mapping[!(pathway_mapping$V2 %in% child_less_parents_filtered), ]
    }
    
  }
  return(pathway_mapping)
}

# Function for getting mama and parents
get_filtered_pathway_names <- function(relation_table, starting_id){
  # get all of the parents of the starting ID
  all_parents <- get_parents(relation_table, starting_id)
  return(all_parents)
}
get_parents <- function(relation_table, starting_id){
  # get all of the parents of the starting ID
  parents <- as.character(relation_table[!is.na(relation_table$V2) & !is.na(relation_table$V1) & relation_table$V2 == starting_id, 'V1'])
  # these parents are all family
  family <- parents
  # see if there were any parents
  if(length(parents) > 0){
    # if there were parents, we need to get their parents as well
    for(parent in parents){
      # get the grandparents and add these to the family
      grand_parents <- get_parents(relation_table, parent)
      family <- c(family, grand_parents)
    }
  }
  return(family)
}

pathways_to_trees <- function(relation_table){
  # first get the biggest parents, the terms which don't are not children
  super_parents <- child_less_parents <- setdiff(as.character(relation_table$V1), as.character(relation_table$V2))
  # we must do a tree per super parent
  super_parent_tree_list <- list()
  # put in the work for each super parent
  for(super_parent in super_parents){
    # fetch complete list with attribute
    super_parent_node <- get_children_lists(relation_table, super_parent, 10)
    # set attributes for super parent itself
    class(super_parent_node) <- 'dendrogram'
    # add to list of superparents
    super_parent_tree_list[[super_parent]] <- super_parent_node
  }
  return(super_parent_tree_list)
}

get_children_lists <- function(relation_table, term, height){
  # get the children of the term
  child_rows <- relation_table[as.character(relation_table$V1) == as.character(term), ]
  # set up the node
  node <- list()
  # if there are no children, this is a leaf
  if(nrow(child_rows) == 0){
    # set up as leaf
    attributes(node) <- list(members=1, h=0, edgetext=term, label=term, leaf=T)
  }
  else{
    children <- child_rows$V2
    # put in a list the results
    i <- 1
    for(child in children){
      # grab the child lists
      child_node <- get_children_lists(relation_table, child, height-1)
      # add this to the node
      node[[i]] <- child_node
      i <- i + 1
    }
    # add attributes to non-leaf node
    attributes(node) <- list(members=length(children),height=height,edgetext=term)
  }
  return(node)
}

add_pathway_levels <- function(relation_table){
  # first get the biggest parents, the terms which don't are not children
  super_parents <- child_less_parents <- setdiff(as.character(relation_table$V1), as.character(relation_table$V2))
  # we must do a tree per super parent
  super_parent_tree_list <- list()
  # put in the work for each super parent
  for(super_parent in super_parents){
    # fetch complete list with attribute
    super_parent_node <- get_children_pathway_levels(relation_table, super_parent, 0)
    # add to list of superparents
    super_parent_tree_list[[super_parent]] <- super_parent_node
  }
  return(super_parent_tree_list)
}

get_children_pathway_levels <- function(relation_table, term, level){
  # get the children of the term
  child_rows <- relation_table[as.character(relation_table$V1) == as.character(term), ]
  # set up the 
  current_level <- data.frame(term=c(term), level=c(level))
  # if there are no children, this is a leaf
  if(nrow(child_rows) > 0){
    children <- child_rows$V2
    for(child in children){
      # grab the child lists
      child_node <- get_children_pathway_levels(relation_table, child, level+1)
      # add to our level
      current_level <- rbind(current_level, child_node)
    }
  }
  return(current_level)
}

plot_de_gene_uniqueness_condition <- function(base_mast_output_path, marked_singles=T, condition_combinations=c('UTBaseline', 'UTt24h', 'UTt8w', 'Baselinet24h', 'Baselinet8w', 't24ht8w'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # we want to get all the DE and where they were significant
  gene_conditions_df <- NULL
  # so, check each condition
  for(condition in condition_combinations){
    # we'll store the genes
    genes_condition <- c()
    # then check each cell type
    for(cell_type in cell_types_to_use){
      # paste together the filename
      full_mast_path <- paste(base_mast_output_path, cell_type, condition, '.tsv', sep = '')
      # get the genes
      genes_cell_type <- get_de_genes_from_mast_file(full_mast_path, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping)
      # add to the list of genes for this condition
      genes_condition <- c(genes_condition, genes_cell_type)
    }
    # make the genes unique, as the cell types have overlap
    genes_condition <- unique(genes_condition)
    # set as a dataframe
    condition_df <- data.frame(genes_condition, rep(T, times=length(genes_condition)), stringsAsFactors = F)
    # set the condition as the column name
    colnames(condition_df) <- c('gene', condition)
    # merge with existing gene and conditions
    if(is.null(gene_conditions_df)){
      gene_conditions_df <- condition_df
    }
    else{
      gene_conditions_df <- merge(gene_conditions_df, condition_df, by = 'gene', all = T)
    }
  }
  # genes that were not found before merging, will be NA, those are not found and thus of course F
  gene_conditions_df[is.na(gene_conditions_df)] <- F
  # disregard the gene column (makes the applies easier)
  gene_conditions_df <- gene_conditions_df[, condition_combinations]
  # get for each genes how many conditions it was significant in
  number_of_times_sig <- apply(gene_conditions_df, 1, sum)
  # turn into plot df
  plot_df <- NULL
  # we need a bar for unique number that the DE genes were unique
  for(number_sig in unique(number_of_times_sig)){
    # get the number of times this was the case
    number_times_this_sig <- sum(number_of_times_sig == number_sig)
    # create appropriate row
    row <- data.frame(number_sig=number_sig, number_times_this_sig=number_times_this_sig, stringsAsFactors = F)
    # add or set
    if(is.null(plot_df)){
      plot_df <- row
    }
    else{
      plot_df <- rbind(plot_df, row)
    }
  }
  # set the condition of the numbers
  plot_df$condition <- 'mixed'
  # for the singles, so in one condition only, we might want to see the proportions
  if(marked_singles){
    # remove the singles, as we're overwriting those
    plot_df <- plot_df[plot_df$number_sig != 1, ]
    # subset to the singles
    gene_conditions_df_singles <- gene_conditions_df[number_of_times_sig == 1, ]
    # now that we have only the singles, we can sum over the columns, to get the number of genes unique to the condition
    number_unique_per_condition <- apply(gene_conditions_df_singles, 2, sum)
    # make that into a df
    plot_df_uniques <- data.frame(number_sig=rep(1, times=length(condition_combinations)), number_times_this_sig=number_unique_per_condition, condition=condition_combinations)
    # add to the current plot
    plot_df <- rbind(plot_df, plot_df_uniques)
  }
  # set the order I like for the legend, but setting the factor order
  plot_df$condition <- factor(plot_df$condition, levels=c('mixed', condition_combinations))
  # grab the colours
  cc <- get_color_coding_dict()
  # add the 'mixed' condition
  cc[['mixed']] <- 'gray'
  fillScale <- scale_fill_manual(name = "condition",values = unlist(cc[c(condition_combinations, 'mixed')]))
  # make the plot finally
  p <- ggplot(plot_df, aes(fill=condition, y=number_times_this_sig, x=number_sig)) +
    geom_bar(position='stack', stat='identity') +
    labs(x='number of conditions a gene is differentially expressed in', y='Number of significant DE genes') +
    ggtitle('overlap of DE genes in condition combinations') +
    labs(fill = "Found in") +
    fillScale

  return(p)
}

plot_de_gene_uniqueness_celltype <- function(base_mast_output_path, marked_singles=T, condition_combinations=c('UTBaseline', 'UTt24h', 'UTt8w', 'Baselinet24h', 'Baselinet8w', 't24ht8w'), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # we want to get all the DE and where they were significant
  gene_cell_type_df <- NULL
  # so, check each condition
  for(cell_type in cell_types_to_use){
    # we'll store the genes
    genes_cell_type <- c()
    # then check each cell type
    for(condition in condition_combinations){
      # paste together the filename
      full_mast_path <- paste(base_mast_output_path, cell_type, condition, '.tsv', sep = '')
      # get the genes
      genes_condition <- get_de_genes_from_mast_file(full_mast_path, pval_column=pval_column, sig_pval=sig_pval, max=max, max_by_pval=max_by_pval, only_positive=only_positive, only_negative=only_negative, lfc_column=lfc_column, to_ens=to_ens, symbols.to.ensg.mapping=symbols.to.ensg.mapping)
      # add to the list of genes for this condition
      genes_cell_type <- c(genes_cell_type, genes_condition)
    }
    # make the genes unique, as the cell types have overlap
    genes_cell_type <- unique(genes_cell_type)
    # set as a dataframe
    cell_type_df <- data.frame(genes_cell_type, rep(T, times=length(genes_cell_type)), stringsAsFactors = F)
    # set the condition as the column name
    colnames(cell_type_df) <- c('gene', cell_type)
    # merge with existing gene and conditions
    if(is.null(gene_cell_type_df)){
      gene_cell_type_df <- cell_type_df
    }
    else{
      gene_cell_type_df <- merge(gene_cell_type_df, cell_type_df, by = 'gene', all = T)
    }
  }
  # genes that were not found before merging, will be NA, those are not found and thus of course F
  gene_cell_type_df[is.na(gene_cell_type_df)] <- F
  # disregard the gene column (makes the applies easier)
  gene_cell_type_df <- gene_cell_type_df[, cell_types_to_use]
  # get for each genes how many conditions it was significant in
  number_of_times_sig <- apply(gene_cell_type_df, 1, sum)
  # turn into plot df
  plot_df <- NULL
  # we need a bar for unique number that the DE genes were unique
  for(number_sig in unique(number_of_times_sig)){
    # get the number of times this was the case
    number_times_this_sig <- sum(number_of_times_sig == number_sig)
    # create appropriate row
    row <- data.frame(number_sig=number_sig, number_times_this_sig=number_times_this_sig, stringsAsFactors = F)
    # add or set
    if(is.null(plot_df)){
      plot_df <- row
    }
    else{
      plot_df <- rbind(plot_df, row)
    }
  }
  # set the condition of the numbers
  plot_df$cell_type <- 'mixed'
  # for the singles, so in one condition only, we might want to see the proportions
  if(marked_singles){
    # remove the singles, as we're overwriting those
    plot_df <- plot_df[plot_df$number_sig != 1, ]
    # subset to the singles
    gene_cell_type_df_singles <- gene_cell_type_df[number_of_times_sig == 1, ]
    # now that we have only the singles, we can sum over the columns, to get the number of genes unique to the condition
    number_unique_per_cell_type <- apply(gene_cell_type_df_singles, 2, sum)
    # make that into a df
    plot_df_uniques <- data.frame(number_sig=rep(1, times=length(cell_types_to_use)), number_times_this_sig=number_unique_per_cell_type, cell_type=cell_types_to_use)
    # add to the current plot
    plot_df <- rbind(plot_df, plot_df_uniques)
  }
  # set the order I like for the legend, but setting the factor order
  plot_df$cell_type <- factor(plot_df$cell_type, levels=c('mixed', cell_types_to_use))
  # grab the colours
  cc <- get_color_coding_dict()
  # add the 'mixed' condition
  cc[['mixed']] <- 'gray'
  fillScale <- scale_fill_manual(name = "cell_type",values = unlist(cc[c(cell_types_to_use, 'mixed')]))
  # make the plot finally
  p <- ggplot(plot_df, aes(fill=cell_type, y=number_times_this_sig, x=number_sig)) +
    geom_bar(position='stack', stat='identity') +
    labs(x='number of cell types a gene is differentially expressed in', y='Number of significant DE genes') +
    ggtitle('overlap of DE genes in cell types') +
    labs(fill = "Found in") +
    fillScale
  
  return(p)
}


get_de_genes_from_mast_file <- function(mast_full_file_path, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  significant_genes <- c()
  try({
    # read the mast output
    mast <- read.table(mast_full_file_path, header=T, row.names = 1)
    # filter to only include the significant results
    mast <- mast[mast[[pval_column]] <= 0.05, ]
    # filter for only the positive lfc if required
    if(only_positive){
      mast <- mast[mast[[lfc_column]] < 0, ]
    }
    # filter for only the positive lfc if required
    if(only_negative){
      mast <- mast[mast[[lfc_column]] > 0, ]
    }
    # confine in some way if reporting a max number of genes
    if(!is.null(max)){
      # by p if required
      if(max_by_pval){
        mast <- mast[order(mast[[pval_column]]), ]
      }
      # by lfc otherwise
      else{
        mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        # if we only have the upregulated genes, order the other way around
        if(only_positive){
          mast <- mast[order(mast[[lfc_column]], decreasing = F), ]
        }
      }
      # subset to the number we requested if max was set
      mast <- mast[1:max,]
    }
    # grab the genes from the column names
    genes <- rownames(mast)
    # convert the symbols to ensemble IDs
    if (to_ens) {
      mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
      mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
      genes <- mapping[match(genes, mapping$V2),"V1"]
    }
    # otherwise change the Seurat replacement back
    else{
      #genes <- gsub("-", "_", genes)
    }
    significant_genes <- genes
  })
  return(significant_genes)
}

get_average_gene_expression_per_ct_and_tp <- function(seurat_object, condition.column = 'timepoint.final', cell.type.column = 'cell_type_lowerres', cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC"), conditions=c('UT', 'Baseline', 't24h', 't8w'), assay='RNA'){
  exp_df <- NULL
  # calculate for each condition
  for(condition in conditions){
    # subset to just the cells of this condition
    seurat_object_condition <- seurat_object[,seurat_object@meta.data[condition.column] == condition]
    # calculate for each cell_type
    for(cell_type in cell_types_to_use){
      # subset to just the cells of the cell type
      seurat_object_cell_type <- seurat_object_condition[,seurat_object_condition@meta.data[cell.type.column] == cell_type]
      # calculate the relevant matrix from the relevant assay
      exp_df_ct_cond <- NULL
      if(assay == 'RNA'){
        DefaultAssay(seurat_object_cell_type) <- 'RNA'
        averages <- apply(seurat_object_cell_type$RNA@data, 1, mean)
        exp_df_ct_cond <- data.frame(condition=rep(condition, times = length(averages)), cell_type=rep(cell_type, times = length(averages)), gene=rownames(seurat_object_cell_type$RNA@data), average=averages)
      }
      else if(assay == 'SCT'){
        DefaultAssay(seurat_object_cell_type) <- 'SCT'
        averages <- apply(seurat_object_cell_type$SCT@counts, 1, mean)
        exp_df_ct_cond <- data.frame(condition=rep(condition, times = length(averages)), cell_type=rep(cell_type, times = length(averages)), gene=rownames(seurat_object_cell_type$SCT@counts), average=averages)
      }
      # paste to the overall df
      if(is.null(exp_df)){
        exp_df <- exp_df_ct_cond
      }
      else{
        exp_df <- rbind(exp_df, exp_df_ct_cond)
      }
    }
  }
  return(exp_df)
}

avg_exp_table_to_hm_table <- function(expression_table){
  # initialise table
  hm_table <- NULL
  # go through the conditions
  for(condition in unique(expression_table$condition)){
    # get the expression for that table
    expression_table_cond <- expression_table[expression_table$condition == condition, c('gene', 'average')]
    # set colnames so that we can merge these later
    colnames(expression_table_cond) <- c('gene', condition)
    # convert to data.table for efficient merging
    expression_table_cond <- data.table(expression_table_cond)
    # try to merge if necessary
    if(is.null(hm_table)){
      hm_table <- expression_table_cond
    }
    else{
      hm_table <- merge(hm_table, expression_table_cond, by='gene')
    }
  }
  # convert back to regular dataframe
  hm_table <- data.frame(hm_table)
  # set rownames
  rownames(hm_table) <- hm_table$gene
  # remove the old gene column
  hm_table$gene <- NULL
  return(hm_table)
}

####################
# Main Code        #
####################

# get the locations of the DE output
#mast_output_prepend <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/stemi_v'
mast_output_prepend <- '/data/cardiology/differential_expression/MAST/results/stemi_v'
#mast_output_append <- '_paired_lores_20200707/rna/'
mast_output_append <- '_paired_lores_lfc025minpct01ncountrna_20210301/rna/'
mast_output_append_lfc01 <- '_paired_lores_lfc01minpct01ncountrna_20210301/rna/'
# write the location of the combined output
#mast_meta_output_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_20200707/rna/'
#mast_meta_output_loc <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_20200707/rna/'
#mast_meta_output_loc <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc025minpct01ncountrna_20200707/rna/'
mast_meta_output_loc_lfc01 <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc01minpct01ncountrna_20210301/rna/'

# write meta output
write_meta_mast(mast_output_prepend, mast_output_append, mast_meta_output_loc)
write_meta_mast(mast_output_prepend, mast_output_append_lfc01, mast_meta_output_loc_lfc01)

 # check the overlap between chemistries
overlap_venn_loc <- '/data/cardiology/differential_expression/MAST/overlap/stemi_meta_paired_lores_lfc01minpct01ncountrna_20201209/rna/'
write_version_chemistry_overlap(mast_meta_output_loc_lfc01, overlap_venn_loc)

# mapping of gene symbol to ensemble id
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
# set the location to write the significant genes
sig_output_loc <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc025minpct01_20200707_ensid/rna/'
sig_output_loc_gs <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc025minpct01_20200707/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc, sig_output_loc, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc, sig_output_loc_gs, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc025minpct01_20200707_ensid/rna/'
sig_up_output_loc_gs <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc025minpct01_20200707/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_up_output_loc, only_positive = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc, sig_up_output_loc_gs, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc025minpct01_20200707_ensid/rna/'
sig_down_output_loc_gs <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc025minpct01_20200707/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_down_output_loc, only_negative = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc, sig_down_output_loc_gs, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# for the lfc01 output as well
sig_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20210301_ensid/rna/'
sig_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20210301/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_output_loc_lfc01, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_output_loc_gs_lfc01, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20210301_ensid/rna/'
sig_up_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20210301/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_up_output_loc_lfc01, only_positive = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_up_output_loc_gs_lfc01, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20210301_ensid/rna/'
sig_down_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20210301/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_down_output_loc_lfc01, only_negative = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_down_output_loc_gs_lfc01, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)

# set the location for the significant genes that were upregulated
sig_up_output_loc_lfc01_v2 <- '/data/cardiology/differential_expression/sigs_pos/v2_paired_lores_lfc01minpct01_20210301_ensid/rna/'
sig_up_output_loc_gs_lfc01_v2 <- '/data/cardiology/differential_expression/sigs_pos/v2_paired_lores_lfc01minpct01_20210301/rna/'
# write the significantly upregulated genes
get_significant_genes(paste(mast_output_prepend, '2', mast_output_append_lfc01, sep = ''), sig_up_output_loc_lfc01_v2, only_positive = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping, pval_column='p_val_adj', lfc_column = 'avg_logFC')
get_significant_genes(paste(mast_output_prepend, '2', mast_output_append_lfc01, sep = ''), sig_up_output_loc_gs_lfc01_v2, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping, pval_column='p_val_adj', lfc_column = 'avg_logFC')


# get the location of the pathways
#pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/meta_paired_lores_lfc01minpct01_20200707_ensid_all/rna/'
pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/meta_paired_lores_lfc01minpct01_20201209_ensid_all/rna/'
# write the combined pathway file
pathway_up_df <- get_pathway_table(pathway_up_output_loc, append = '_sig_up_pathways.txt')
pathway_up_df[pathway_up_df==0] <- 350
# get the df limited by top pathways of upregulated genes
pathway_up_df_top_3 <- get_top_pathways(pathway_up_df, 3, T)
pathway_up_df_top_5 <- get_top_pathways(pathway_up_df, 5, T)

# show pathways
cc <- get_color_coding_dict()
colors_cond <- rep(c(cc[['UTBaseline']],cc[['UTt24h']],cc[['UTt8w']],cc[['Baselinet24h']],cc[['Baselinet8w']],cc[['t24ht8w']]), times = 6)
colors_ct <- c(rep(cc[['B']], times=6),rep(cc[['CD4T']], times=6),rep(cc[['CD8T']], times=6),rep(cc[['DC']], times=6),rep(cc[['monocyte']], times=6),rep(cc[['NK']], times=6))
colors_m <- cbind(colors_ct, colors_cond)
colnames(colors_m) <- c('celltype',
                        'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all.png', width = 1600, height = 1200)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all.pdf', width = 16, height = 12)
heatmap.3(t(as.matrix(pathway_up_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))
dev.off()

# get the location of the pathways
v2_pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/v2_paired_lores_lfc01minpct01_20201209_ensid/rna/'
# write the combined pathway file
v2_pathway_up_df <- get_pathway_table(v2_pathway_up_output_loc, append = '_sig_up_pathways.txt')
v2_pathway_up_df[v2_pathway_up_df==0] <- 350
# get the df limited by top pathways of upregulated genes
v2_pathway_up_df_top_3 <- get_top_pathways(v2_pathway_up_df, 3, T)
#png('/data/cardiology/plots/pathways/v2pathway_ensid_up_all.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2pathway_ensid_up_all.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))
dev.off()

# get the location of the pathways
pathway_down_output_loc <- '/data/cardiology/pathways/sigs_neg/meta_paired_lores_lfc01minpct01_20201209_ensid_all/rna/'
# write the combined pathway file
pathway_down_df <- get_pathway_table(pathway_down_output_loc, append = '_sig_down_pathways.txt')
pathway_down_df[pathway_down_df==0] <- 500
# get the df limited by top pathways of upregulated genes
pathway_down_df_top_3 <- get_top_pathways(pathway_down_df, 3, T)

heatmap.3(t(as.matrix(pathway_down_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))
dev.off()


# get the location of the pathways
v2_pathway_down_output_loc <- '/data/cardiology/pathways/sigs_neg/v2_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
v2_pathway_down_df <- get_pathway_table(v2_pathway_down_output_loc, append = '_signeg.txt')
v2_pathway_down_df[v2_pathway_down_df==0] <- 400
# get the df limited by top pathways of upregulated genes
v2_pathway_down_df_top_3 <- get_top_pathways(v2_pathway_down_df, 3, T)

heatmap.3(t(as.matrix(v2_pathway_down_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# create subset heatmaps of conditions
pathway_up_df_baselinet24h <- pathway_up_df[, colnames(pathway_up_df)[grep('Baselinet24h', colnames(pathway_up_df))]]
pathway_up_df_baselinet24h_top_10 <- get_top_pathways(pathway_up_df_baselinet24h, 10, T)
cc <- get_color_coding_dict()
t24h_vary_colors_cond <- rep(c(cc[['Baselinet24h']]), times = 6)
t24h_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t24h_vary_colors_m <- cbind(t24h_vary_colors_ct, t24h_vary_colors_cond)
colnames(t24h_vary_colors_m) <- c('celltype',
                                  'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_baselinet24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_baselinet24h_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_baselinet24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t0 vs t24h')
dev.off()

# now down
pathway_down_df_baselinet24h <- pathway_down_df[, colnames(pathway_down_df)[grep('Baselinet24h', colnames(pathway_down_df))]]
pathway_down_df_baselinet24h_top_10 <- get_top_pathways(pathway_down_df_baselinet24h, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_baselinet24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_baselinet24h_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_baselinet24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(29,11), main = 'downregulated pathways t0 vs t24h')
dev.off()
# baseline vs t8w
pathway_up_df_baselinet8w <- pathway_up_df[, colnames(pathway_up_df)[grep('Baselinet8w', colnames(pathway_up_df))]]
pathway_up_df_baselinet8w_top_10 <- get_top_pathways(pathway_up_df_baselinet8w, 10, T)
baselinet8w_vary_colors_cond <- rep(c(cc[['Baselinet8w']]), times = 6)
baselinet8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
baselinet8w_vary_colors_m <- cbind(baselinet8w_vary_colors_ct, baselinet8w_vary_colors_cond)
colnames(baselinet8w_vary_colors_m) <- c('celltype',
                                  'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_baselinet8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_baselinet8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_baselinet8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(baselinet8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t0 vs t8w')
dev.off()
# now down
pathway_down_df_baselinet8w <- pathway_down_df[, colnames(pathway_down_df)[grep('Baselinet8w', colnames(pathway_down_df))]]
pathway_down_df_baselinet8w_top_10 <- get_top_pathways(pathway_down_df_baselinet8w, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_baselinet8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_baselinet8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_baselinet8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(baselinet8w_vary_colors_m), margins=c(29,11), main = 'downregulated pathways t0 vs t8w')
dev.off()
# t24h vs t8w
pathway_up_df_t24ht8w <- pathway_up_df[, colnames(pathway_up_df)[grep('t24ht8w', colnames(pathway_up_df))]]
pathway_up_df_t24ht8w_top_10 <- get_top_pathways(pathway_up_df_t24ht8w, 10, T)
t24ht8w_vary_colors_cond <- rep(c(cc[['t24ht8w']]), times = 6)
t24ht8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t24ht8w_vary_colors_m <- cbind(t24ht8w_vary_colors_ct, t24ht8w_vary_colors_cond)
colnames(t24ht8w_vary_colors_m) <- c('celltype',
                                  'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_t24ht8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_t24ht8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_t24ht8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t24h vs t8w')
dev.off()
# now down
pathway_down_df_t24ht8w <- pathway_down_df[, colnames(pathway_down_df)[grep('t24ht8w', colnames(pathway_down_df))]]
pathway_down_df_t24ht8w_top_10 <- get_top_pathways(pathway_down_df_t24ht8w, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_t24ht8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_t24ht8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_t24ht8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(29,11), main = 'downregulated pathways t24h vs t8w')
dev.off()

# UT vs t8w
pathway_up_df_utt8w <- pathway_up_df[, colnames(pathway_up_df)[grep('UTt8w', colnames(pathway_up_df))]]
pathway_up_df_utt8w_top_10 <- get_top_pathways(pathway_up_df_utt8w, 10, T)
utt8w_vary_colors_cond <- rep(c(cc[['UTt8w']]), times = 6)
utt8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
utt8w_vary_colors_m <- cbind(utt8w_vary_colors_ct, utt8w_vary_colors_cond)
colnames(utt8w_vary_colors_m) <- c('celltype',
                                     'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_UTt8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_UTt8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_utt8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t8w')
dev.off()

# UT vs Baseline
pathway_up_df_utbaseline <- pathway_up_df[, colnames(pathway_up_df)[grep('UTBaseline', colnames(pathway_up_df))]]
pathway_up_df_utbaseline_top_10 <- get_top_pathways(pathway_up_df_utbaseline, 10, T)
utbaseline_vary_colors_cond <- rep(c(cc[['UTBaseline']]), times = 6)
utbaseline_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
utbaseline_vary_colors_m <- cbind(utbaseline_vary_colors_ct, utbaseline_vary_colors_cond)
colnames(utbaseline_vary_colors_m) <- c('celltype',
                                   'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_UTBaseline_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_UTBaseline_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_utbaseline_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utbaseline_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t0')
dev.off()

# UT vs t24h
pathway_up_df_utt24h <- pathway_up_df[, colnames(pathway_up_df)[grep('UTt24h', colnames(pathway_up_df))]]
pathway_up_df_utt24h_top_10 <- get_top_pathways(pathway_up_df_utt24h, 10, T)
utt24h_vary_colors_cond <- rep(c(cc[['UTt24h']]), times = 6)
utt24h_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
utt24h_vary_colors_m <- cbind(utt24h_vary_colors_ct, utt24h_vary_colors_cond)
colnames(utt24h_vary_colors_m) <- c('celltype',
                                        'condition')
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_ut24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_up_all_ut24h_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_up_df_utt24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt24h_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t24h')
dev.off()

# UT downs
pathway_down_df_utt8w <- pathway_down_df[, colnames(pathway_down_df)[grep('UTt8w', colnames(pathway_down_df))]]
max(pathway_down_df_utt8w)
pathway_down_df_utt8w[pathway_down_df_utt8w == 0] <- 400
pathway_down_df_utt8w_top_10 <- get_top_pathways(pathway_down_df_utt8w, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_utt8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_utt8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_utt8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt8w_vary_colors_m), margins=c(29,11), main = 'downregulated pathways HC vs t8w')
dev.off()

pathway_down_df_utt24h <- pathway_down_df[, colnames(pathway_down_df)[grep('UTt24h', colnames(pathway_down_df))]]
max(pathway_down_df_utt24h)
pathway_down_df_utt24h[pathway_down_df_utt24h == 0] <- 300
pathway_down_df_utt24h_top_10 <- get_top_pathways(pathway_down_df_utt24h, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_ut24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_ut24h_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_utt24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt24h_vary_colors_m), margins=c(29,11), main = 'downregulated pathways HC vs t24h')
dev.off()

pathway_down_df_utbaseline <- pathway_down_df[, colnames(pathway_down_df)[grep('UTBaseline', colnames(pathway_down_df))]]
max(pathway_down_df_utbaseline)
pathway_down_df_utbaseline[pathway_down_df_utbaseline == 0] <- 400
pathway_down_df_utbaseline_top_10 <- get_top_pathways(pathway_down_df_utbaseline, 10, T)
#png('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_utbaseline_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/meta_pathway_ensid_down_all_utbaseline_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(pathway_down_df_utbaseline_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utbaseline_vary_colors_m), margins=c(29,11), main = 'downregulated pathways HC vs t0')
dev.off()

# v2 specific
v2_pathway_up_df_baselinet24h <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('Baselinet24h', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_baselinet24h_top_10 <- get_top_pathways(v2_pathway_up_df_baselinet24h, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_baselinet24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_baselinet24h_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_baselinet24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t0 vs t24h')
dev.off()
# baseline vs t8w
v2_pathway_up_df_baselinet8w <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('Baselinet8w', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_baselinet8w_top_10 <- get_top_pathways(v2_pathway_up_df_baselinet8w, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_baselinet8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_baselinet8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_baselinet8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(baselinet8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t0 vs t8w')
dev.off()
# t24h vs t8w
v2_pathway_up_df_t24ht8w <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('t24ht8w', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_t24ht8w_top_10 <- get_top_pathways(v2_pathway_up_df_t24ht8w, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_t24ht8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_t24ht8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_t24ht8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways t24h vs t8w')
dev.off()
# UT vs t8w
v2_pathway_up_df_utt8w <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('UTt8w', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_utt8w_top_10 <- get_top_pathways(v2_pathway_up_df_utt8w, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_UTt8w_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_UTt8w_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_utt8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt8w_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t8w')
dev.off()
# UT vs Baseline
v2_pathway_up_df_utbaseline <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('UTBaseline', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_utbaseline_top_10 <- get_top_pathways(v2_pathway_up_df_utbaseline, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_UTBaseline_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_UTBaseline_top_10.pdf', width = 12, height = 9)
heatmap.3(t(as.matrix(v2_pathway_up_df_utbaseline_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utbaseline_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t0')
dev.off()
# UT vs t24h
v2_pathway_up_df_utt24h <- v2_pathway_up_df[, colnames(v2_pathway_up_df)[grep('UTt24h', colnames(v2_pathway_up_df))]]
v2_pathway_up_df_utt24h_top_10 <- get_top_pathways(v2_pathway_up_df_utt24h, 10, T)
#png('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_ut24h_top_10.png', width = 1200, height = 900)
pdf('/data/cardiology/plots/pathways/v2_pathway_ensid_up_all_ut24h_top_10.pdf', width=12, height=9)
heatmap.3(t(as.matrix(v2_pathway_up_df_utt24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt24h_vary_colors_m), margins=c(29,11), main = 'upregulated pathways HC vs t24h')
dev.off()


# create heatmap of Baseline vs t24h
cc <- get_color_coding_dict()
stemi_de_table <- get_combined_meta_de_table(mast_meta_output_loc_lfc01)
# do per combination
stemi_de_table_baseline_t24h <- stemi_de_table[, colnames(stemi_de_table)[grep('Baselinet24h', colnames(stemi_de_table))]]
stemi_de_table_baseline_t24h_vary <- get_top_vary_genes(stemi_de_table_baseline_t24h, use_ct = F, use_tp = T, top_so_many = 50, use_dynamic_sd = T, timepoints=c("Baselinet24h"))
deg_meta_stemi_de_table_baseline_t24h_ct_vary <- stemi_de_table_baseline_t24h[(rownames(stemi_de_table_baseline_t24h) %in% stemi_de_table_baseline_t24h_vary), ]
t24h_vary_colors_cond <- rep(c(cc[['Baselinet24h']]), times = 6)
t24h_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t24h_vary_colors_m <- cbind(t24h_vary_colors_ct, t24h_vary_colors_cond)
colnames(t24h_vary_colors_m) <- c('celltype',
                        'condition')
heatmap.3(t(as.matrix(deg_meta_stemi_de_table_baseline_t24h_ct_vary)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(5,11), main = 't0/t24h DE lfc varied')
heatmap.3(t(as.matrix(stemi_de_table_baseline_t24h)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(5,11), labCol=NA, main = 't0/t24h DE lfc')

# create heatmp of Baseline vs t8w
stemi_de_table_baseline_t8w <- stemi_de_table[, colnames(stemi_de_table)[grep('Baselinet8w', colnames(stemi_de_table))]]
stemi_de_table_baseline_t8w_vary <- get_top_vary_genes(stemi_de_table_baseline_t8w, use_ct = F, use_tp = T, top_so_many = 50, use_dynamic_sd = T, timepoints=c("Baselinet8w"))
deg_meta_stemi_de_table_baseline_t8w_ct_vary <- stemi_de_table_baseline_t8w[(rownames(stemi_de_table_baseline_t8w) %in% stemi_de_table_baseline_t8w_vary), ]
t8w_vary_colors_cond <- rep(c(cc[['Baselinet8w']]), times = 6)
t8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t8w_vary_colors_m <- cbind(t8w_vary_colors_ct, t8w_vary_colors_cond)
colnames(t8w_vary_colors_m) <- c('celltype',
                                  'condition')
heatmap.3(t(as.matrix(deg_meta_stemi_de_table_baseline_t8w_ct_vary)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t8w_vary_colors_m), margins=c(5,11), main = 't0/t8w DE lfc varied')
heatmap.3(t(as.matrix(stemi_de_table_baseline_t8w)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t8w_vary_colors_m), margins=c(5,11), labCol=NA, main = 't0/t8w DE lfc')

# create heatmp of t24h vs t8w
stemi_de_table_t24h_t8w <- stemi_de_table[, colnames(stemi_de_table)[grep('t24ht8w', colnames(stemi_de_table))]]
stemi_de_table_t24h_t8w_vary <- get_top_vary_genes(stemi_de_table_t24h_t8w, use_ct = F, use_tp = T, top_so_many = 50, use_dynamic_sd = T, timepoints=c("t24ht8w"))
deg_meta_stemi_de_table_t24h_t8w_ct_vary <- stemi_de_table_t24h_t8w[(rownames(stemi_de_table_t24h_t8w) %in% stemi_de_table_t24h_t8w_vary), ]
t24ht8w_vary_colors_cond <- rep(c(cc[['t24ht8w']]), times = 6)
t24ht8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t24ht8w_vary_colors_m <- cbind(t24ht8w_vary_colors_ct, t24ht8w_vary_colors_cond)
colnames(t24ht8w_vary_colors_m) <- c('celltype',
                                 'condition')
heatmap.3(t(as.matrix(deg_meta_stemi_de_table_t24h_t8w_ct_vary)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(5,11), main = 't24h/t8w DE lfc varied')
heatmap.3(t(as.matrix(stemi_de_table_t24h_t8w)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(5,11), labCol=NA, main = 't24h/t8w DE lfc')

# create the full table and heatmap
stemi_de_table_andut <- get_combined_meta_de_table(mast_meta_output_loc_lfc01, timepoints=c("UTBaseline", "UTt24h", "UTt8w", "Baselinet24h", "Baselinet8w", "t24ht8w"))
stemi_de_table_andut_vary <- get_top_vary_genes(stemi_de_table_andut, use_ct = F, use_tp = F, top_so_many = 50, use_dynamic_sd = T, timepoints=c("UTBaseline", "UTt24h", "UTt8w", "Baselinet24h", "Baselinet8w", "t24ht8w"))
deg_meta_stemi_de_table_andut_vary <- stemi_de_table_andut[(rownames(stemi_de_table_andut) %in% stemi_de_table_andut_vary), ]
colors_cond <- rep(c(cc[['UTBaseline']],cc[['UTt24h']],cc[['UTt8w']],cc[['Baselinet24h']],cc[['Baselinet8w']],cc[['t24ht8w']]), times = 6)
colors_ct <- c(rep(cc[['B']], times=6),rep(cc[['CD4T']], times=6),rep(cc[['CD8T']], times=6),rep(cc[['DC']], times=6),rep(cc[['monocyte']], times=6),rep(cc[['NK']], times=6))
colors_m <- cbind(colors_ct, colors_cond)
colnames(colors_m) <- c('celltype',
                        'condition')
heatmap.3(t(as.matrix(deg_meta_stemi_de_table_andut_vary)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(5,9), main = 'DE lfc varied')
heatmap.3(t(as.matrix(stemi_de_table_andut)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(5,9), labCol=NA, main = 'DE lfc')

# this is the reactome ID for the immune system
immune_system_reactome_id <- 'R-HSA-168256'
# load the pathways
pathways <- read.table('/data/scRNA/pathways/ReactomePathways20210115.tsv', sep='\t', quote='')
# subset to just human to speed up the search
pathways <- pathways[pathways$V3 == 'Homo sapiens', ]
# load the pathway mapping
pathway_mappings <- read.table('/data/scRNA/pathways/ReactomePathwaysRelation.tsv', sep = '\t')
# get the filtered names
filtered_names <- get_filtered_pathway_names(pathways, pathway_mappings, 'R-HSA-168256')
# get the df that is left after filtering
pathway_up_df_filtered <- filter_pathway_df_on_starting_id(pathway_up_df, filtered_names)
# check what is top now
pathway_up_df_filtered_top_5 <- get_top_pathways(pathway_up_df_filtered, 5, T)
heatmap.3(t(as.matrix(pathway_up_df_filtered)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(18,10))



# grab the metadata
meta.data <- read.table('/data/cardiology/metadata/cardio.integrated.20201209.metadata.tsv', sep='\t', header=T, row.names=1)
plot_de_vs_cell_type_numbers(mast_meta_output_loc_lfc01, meta.data, plot_separately = F, proportion = T)
plot_de_vs_cell_type_numbers(mast_meta_output_loc_lfc01, meta.data, plot_separately = F, proportion = F)
plot_de_vs_cell_type_numbers(mast_meta_output_loc_lfc01, meta.data, plot_separately = F, proportion = F, cell_types_to_use = c('monocyte'))

# get a table of the number of DE genes
de_numbers_table <- de_genes_number_to_table(mast_meta_output_loc_lfc01)
de_numbers_table_up <- de_genes_number_to_table(mast_meta_output_loc_lfc01, only_positive = T)
de_numbers_table_down <- de_genes_number_to_table(mast_meta_output_loc_lfc01, only_negative = T)
# plot the number of DE genes
numbers_table_to_plot(de_numbers_table)
numbers_table_to_plot(de_numbers_table, use_groups_dict = F, cols_include = c('UTBaseline', 'UTt24h', 'UTt8w'))
numbers_table_to_plot(de_numbers_table, use_groups_dict = F, cols_include = c('Baselinet24h', 'Baselinet8w', 't24ht8w'), pointless = T)
# only HC vs t0
numbers_table_to_plot(de_numbers_table, use_groups_dict = F, cols_include = c('UTBaseline'))


# specifically for monocytes, check the number of DE genes and the fractional differences of their sub-celltype populations
plot_de_number_vs_subcell_population(mast_meta_output_loc_lfc01, 'monocyte', 'cMono', meta.data, timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', make_absolute=F)
plot_de_number_vs_subcell_population(mast_meta_output_loc_lfc01, 'monocyte', 'cMono', meta.data[meta.data$chem=='V2',], timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', make_absolute=F)
plot_de_number_vs_subcell_population(mast_meta_output_loc_lfc01, 'monocyte', 'cMono', meta.data[meta.data$chem=='V3',], timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', make_absolute=F)

# get the top five for each cell type and cell type and condition in terms of up/down regulated genes
top_fives_de_genes <- get_top_de_genes_per_cond_and_ct(mast_meta_output_loc_lfc01, max_by_pval = F, stims = c('Baseline', 't24h', 't8w'))
# get the monocyte genes
mono_up_or_down_genes <- as.vector(unlist(top_fives_de_genes[, colnames(top_fives_de_genes)[grep('monocyte', colnames(top_fives_de_genes))]]))
# plotting the expression instead of the LFC
v2_exp_loc <- '/data/cardiology/differential_expression/cardio.integrated.v2.20201209.avg.exp.tsv'
v3_exp_loc <- '/data/cardiology/differential_expression/cardio.integrated.v3.20201209.avg.exp.tsv'
# get expression
v2_exp <- read.table(v2_exp_loc, sep = '\t', header = T)
v3_exp <- read.table(v3_exp_loc, sep = '\t', header = T)
# confine to DE genes selected
v2_exp_de <- v2_exp[v2_exp$gene %in% mono_up_or_down_genes,]
v3_exp_de <- v3_exp[v3_exp$gene %in% mono_up_or_down_genes,]
# grab only the monocytes
cell_types_to_use <- c('monocyte')
conditions_to_use <- c('Baseline', 't24h', 't8w')
v2_exp_de <- v2_exp_de[v2_exp_de$cell_type %in% cell_types_to_use, ]
v3_exp_de <- v3_exp_de[v3_exp_de$cell_type %in% cell_types_to_use, ]
v2_exp_de <- v2_exp_de[v2_exp_de$condition %in% conditions_to_use, ]
v3_exp_de <- v3_exp_de[v3_exp_de$condition %in% conditions_to_use, ]
# turn into hm
v2_exp_hm_format <- avg_exp_table_to_hm_table(v2_exp_de)
v3_exp_hm_format <- avg_exp_table_to_hm_table(v3_exp_de)
# replace Baseline with t0
colnames(v2_exp_hm_format) <- gsub('Baseline', 't0', colnames(v2_exp_hm_format))
colnames(v3_exp_hm_format) <- gsub('Baseline', 't0', colnames(v3_exp_hm_format))

# plot
heatmap.3((as.matrix(v2_exp_hm_format)),
          col=rev(brewer.pal(10,"RdBu")), margins=c(6,10), dendrogram = 'none')
heatmap.3((as.matrix(v3_exp_hm_format)),
          col=rev(brewer.pal(10,"RdBu")), margins=c(6,10), dendrogram = 'none')


# plot the sharing of DE genes
plot_DE_sharing_per_celltype('UTBaseline', mast_meta_output_loc_lfc01)
plot_DE_sharing_per_celltype('UTt24h', mast_meta_output_loc_lfc01)
plot_DE_sharing_per_celltype('UTt8w', mast_meta_output_loc_lfc01)
plot_DE_sharing_per_celltype('Baselinet24h', mast_meta_output_loc_lfc01)
plot_DE_sharing_per_celltype('t24ht8w', mast_meta_output_loc_lfc01)
plot_DE_sharing_per_celltype('Baselinet8w', mast_meta_output_loc_lfc01)

plot_DE_sharing_per_celltype('UTBaseline', mast_meta_output_loc_lfc01, only_positive = T)
plot_DE_sharing_per_celltype('UTt24h', mast_meta_output_loc_lfc01, only_positive = T)
plot_DE_sharing_per_celltype('UTt8w', mast_meta_output_loc_lfc01, only_positive = T)
plot_DE_sharing_per_celltype('Baselinet24h', mast_meta_output_loc_lfc01, only_positive = T)
plot_DE_sharing_per_celltype('t24ht8w', mast_meta_output_loc_lfc01, only_positive = T)
plot_DE_sharing_per_celltype('Baselinet8w', mast_meta_output_loc_lfc01, only_positive = T)


mast_lfc01_de_genes <- get_de_genes(mast_output_loc = mast_meta_output_loc_lfc01)
upset(fromList(mast_lfc01_de_genes[['monocyte']]), nsets = length(mast_lfc01_de_genes[['monocyte']]), order.by = 'freq')
mast_lfc01_de_genes_up <- get_de_genes(mast_output_loc = mast_meta_output_loc_lfc01, only_positive = T)
mast_lfc01_de_genes_down <- get_de_genes(mast_output_loc = mast_meta_output_loc_lfc01, only_negative = T)
names(mast_lfc01_de_genes_up[['monocyte']]) <- paste(names(mast_lfc01_de_genes_up[['monocyte']]), '_up')
names(mast_lfc01_de_genes_down[['monocyte']]) <- paste(names(mast_lfc01_de_genes_down[['monocyte']]), '_down')
upset(fromList(append(mast_lfc01_de_genes_up[['monocyte']], mast_lfc01_de_genes_down[['monocyte']])), order.by = 'freq', length(append(mast_lfc01_de_genes_up[['monocyte']], mast_lfc01_de_genes_down[['monocyte']])))

# Load table with pathway
monoUTBase <- read.table('/data/cardiology/pathways/sigs_pos/meta_paired_lores_lfc01minpct01_20201209_ensid_all/rna/monocyteUTBaseline_sig_up_pathways.txt', sep = '\t', header = T, dec = ",")

# Selecting for REACTOME pathways
monoUTBase <- monoUTBase[monoUTBase$Source =='BioSystems: REACTOME', ] 

# Making the variable as numeric
monoUTBase$q.value.FDR.B.H = as.numeric(as.character(monoUTBase$q.value.FDR.B.H))
# match the ID
monoUTBase$matchID <- pathways[match(toupper(as.character(monoUTBase$Name)), toupper(as.character(pathways$V2))), ]$V1
# subset to human pathways
pathway_mappings <- pathway_mappings[startsWith(pathway_mappings$V1, 'R-HSA'), ]

# Fetch MatchIDs
all_mamas_per_starting_id <- list()
for(starting_id in monoUTBase$matchID){
  if(!is.na(starting_id)){
    all_mamas <- get_filtered_pathway_names(pathway_mappings, starting_id)
    all_mamas_per_starting_id[[starting_id]] <- all_mamas
  }
}

# Fetch missing MatchIDs
monoUTBase$matchID[3] <- "R-HSA-168898"
monoUTBase$matchID[4] <- "R-HSA-166058"
"Activated TLR4 signalling" 
"MyD88-independent TLR3/TLR4 cascade"
"TRIF-mediated TLR3/TLR4 signaling"
"RIG-I/MDA5 mediated induction of IFN-alpha/beta pathways"
monoUTBase$matchID[36] <- "R-HSA-6785807"

# get the pathways that are in our dataset or have parents in our dataset
filtered_pathways <- pathway_mapping_filtered_childless(all_mamas_per_starting_id, pathway_mappings)
# filter further to only go from the starting id
#filtered_pathways <- filtered_pathways[filtered_pathways$V1 %in% get_children(pathway_mappings, 'R-HSA-168256') | filtered_pathways$V2 %in% get_children(pathway_mappings, 'R-HSA-168256'), ]

# turn these into treemaps
super_trees <- pathways_to_trees(filtered_pathways)

# find the super parents
super_parent_pathways <- setdiff(as.character(filtered_pathways$V1), as.character(filtered_pathways$V2))

# add for the super parents, another super parent
filtered_pathways <- rbind(data.frame(V1=rep('pathways', times=length(super_parent_pathways)), V2=super_parent_pathways), filtered_pathways)

# get pathway levels
pathway_levels <- add_pathway_levels(filtered_pathways)

# create datatable
my_network <- data.table::data.table(from=as.character(filtered_pathways$V1), to=as.character(filtered_pathways$V2))
#my_familytree <- get_familytree(my_network, "R-HSA-168256")
my_familytree <- get_familytree(my_network, "pathways")
# add the id as title
my_familytree$nodes$title <- my_familytree$nodes$id
# add the actual name as the title
my_familytree$nodes$title <- pathways[match(as.character(my_familytree$nodes$id), as.character(pathways$V1)), ]$V2
# replace NAs
my_familytree$nodes[is.na(my_familytree$nodes$title), ]$title <- 'pathways'
# set the size from the P value
my_familytree$nodes$value <- monoUTBase[match(as.character(my_familytree$nodes$id), as.character(monoUTBase$matchID)), ]$q.value.FDR.B.H
# turn into number
my_familytree$nodes$value <- as.numeric(my_familytree$nodes$value)
# set the empty Ps into '1'
my_familytree$nodes[is.na(my_familytree$nodes$value), ]$value <- 1
# add the P to the title
my_familytree$nodes$title <- paste(my_familytree$nodes$title, ' p=', as.character(my_familytree$nodes$value), sep='')
# reverse size, we want a smaller P to be a bigger value
my_familytree$nodes$value <- 1 - my_familytree$nodes$value
# add the color by group
#my_familytree$nodes$group <- pathway_levels[['R-HSA-168256']][match(my_familytree$nodes$id, pathway_levels[['R-HSA-168256']]$term), ]$level
my_familytree$nodes$group <- pathway_levels[['pathways']][match(my_familytree$nodes$id, pathway_levels[['pathways']]$term), ]$level

# visualize the network
visNetwork(my_familytree$nodes, my_familytree$edges) %>%
  visEdges(arrows = "to")


