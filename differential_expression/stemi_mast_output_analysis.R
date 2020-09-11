####################
# libraries        #
####################

library(RColorBrewer)
library(MetaVolcanoR)
library(ggplot2)
library(data.table)
library(ggpubr)

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
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#E64B50"
  color_coding[["NK"]] <- "#EDBA1B"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
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


plot_de_vs_cell_type_numbers <- function(mast_output_loc, metadata, timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_types_to_use=c("B", "CD4T", "CD8T", "DC", "monocyte", "NK"), pval_column='metap_bonferroni', sig_pval=0.05, plot_separately=T, proportion=F){
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
        geom_point()
    }
  }
  print(table)
  if(!plot_separately){
    # the label is different depending on whether we show the proportion or not
    xlab <- 'nr of cells'
    if(proportion){
      xlab <- 'proportion of cells'
    }
    
    ggplot(table, aes(x=as.numeric(nr_of_cells), y=as.numeric(sig_gene_number), shape=timepoint, color=cell_type)) +
      #scale_x_discrete(breaks = seq(0, 1, by = 0.1)) +
      #scale_y_discrete(breaks = seq(0, 1000, by = 100)) +
      geom_point(size=3) +
      labs(x = xlab, y = 'nr of DE genes', title = 'cells vs nr of DE genes')
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


####################
# Main Code        #
####################

# get the locations of the DE output
#mast_output_prepend <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/stemi_v'
mast_output_prepend <- '/data/cardiology/differential_expression/MAST/results/stemi_v'
#mast_output_append <- '_paired_lores_20200707/rna/'
mast_output_append <- '_paired_lores_lfc025minpct01ncountrna_20200707/rna/'
mast_output_append_lfc01 <- '_paired_lores_lfc01minpct01ncountrna_20200707/rna/'
# write the location of the combined output
#mast_meta_output_loc <- '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_20200707/rna/'
#mast_meta_output_loc <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_20200707/rna/'
mast_meta_output_loc <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc025minpct01ncountrna_20200707/rna/'
mast_meta_output_loc_lfc01 <- '/data/cardiology/differential_expression/MAST/results/stemi_meta_paired_lores_lfc01minpct01ncountrna_20200707/rna/'

# write meta output
write_meta_mast(mast_output_prepend, mast_output_append, mast_meta_output_loc)
write_meta_mast(mast_output_prepend, mast_output_append_lfc01, mast_meta_output_loc_lfc01)

# check the overlap between chemistries
overlap_venn_loc <- '/data/cardiology/differential_expression/MAST/overlap/stemi_meta_paired_lores_lfc01minpct01ncountrna_20200707/rna/'
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
sig_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
sig_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20200707/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_output_loc_lfc01, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_output_loc_gs_lfc01, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
sig_up_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20200707/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_up_output_loc_lfc01, only_positive = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_up_output_loc_gs_lfc01, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc_lfc01 <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
sig_down_output_loc_gs_lfc01 <- '/data/cardiology/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20200707/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc_lfc01, sig_down_output_loc_lfc01, only_negative = T, to_ens = T, symbols.to.ensg.mapping = gene_to_ens_mapping)
get_significant_genes(mast_meta_output_loc_lfc01, sig_down_output_loc_gs_lfc01, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)


# get the location of the pathways
pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
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
heatmap.3(t(as.matrix(pathway_up_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# get the location of the pathways
v2_pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/v2_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
v2_pathway_up_df <- get_pathway_table(v2_pathway_up_output_loc, append = '_sigpos.txt')
v2_pathway_up_df[v2_pathway_up_df==0] <- 350
# get the df limited by top pathways of upregulated genes
v2_pathway_up_df_top_3 <- get_top_pathways(v2_pathway_up_df, 3, T)

heatmap.3(t(as.matrix(v2_pathway_up_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# get the location of the pathways
pathway_down_output_loc <- '/data/cardiology/pathways/sigs_neg/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
pathway_down_df <- get_pathway_table(pathway_down_output_loc, append = '_pathwayneg.txt')
pathway_down_df[pathway_down_df==0] <- 400
# get the df limited by top pathways of upregulated genes
pathway_down_df_top_3 <- get_top_pathways(pathway_down_df, 3, T)

heatmap.3(t(as.matrix(pathway_down_df_top_3)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

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
heatmap.3(t(as.matrix(pathway_up_df_baselinet24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(11,11), main = 'upregulated pathways t0 vs t24h')
# now down
pathway_down_df_baselinet24h <- pathway_down_df[, colnames(pathway_down_df)[grep('Baselinet24h', colnames(pathway_down_df))]]
pathway_down_df_baselinet24h_top_10 <- get_top_pathways(pathway_down_df_baselinet24h, 10, T)
heatmap.3(t(as.matrix(pathway_down_df_baselinet24h_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24h_vary_colors_m), margins=c(11,11), main = 'downregulated pathways t0 vs t24h')
# baseline vs t8w
pathway_up_df_baselinet8w <- pathway_up_df[, colnames(pathway_up_df)[grep('Baselinet8w', colnames(pathway_up_df))]]
pathway_up_df_baselinet8w_top_10 <- get_top_pathways(pathway_up_df_baselinet8w, 10, T)
baselinet8w_vary_colors_cond <- rep(c(cc[['Baselinet8w']]), times = 6)
baselinet8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
baselinet8w_vary_colors_m <- cbind(baselinet8w_vary_colors_ct, baselinet8w_vary_colors_cond)
colnames(baselinet8w_vary_colors_m) <- c('celltype',
                                  'condition')
heatmap.3(t(as.matrix(pathway_up_df_baselinet8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(baselinet8w_vary_colors_m), margins=c(11,11), main = 'upregulated pathways t0 vs t8w')
# now down
pathway_down_df_baselinet8w <- pathway_down_df[, colnames(pathway_down_df)[grep('Baselinet8w', colnames(pathway_down_df))]]
pathway_down_df_baselinet8w_top_10 <- get_top_pathways(pathway_down_df_baselinet8w, 10, T)
heatmap.3(t(as.matrix(pathway_down_df_baselinet8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(baselinet8w_vary_colors_m), margins=c(11,11), main = 'downregulated pathways t0 vs t8w')
# t24h vs t8w
pathway_up_df_t24ht8w <- pathway_up_df[, colnames(pathway_up_df)[grep('t24ht8w', colnames(pathway_up_df))]]
pathway_up_df_t24ht8w_top_10 <- get_top_pathways(pathway_up_df_t24ht8w, 10, T)
t24ht8w_vary_colors_cond <- rep(c(cc[['t24ht8w']]), times = 6)
t24ht8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
t24ht8w_vary_colors_m <- cbind(t24ht8w_vary_colors_ct, t24ht8w_vary_colors_cond)
colnames(t24ht8w_vary_colors_m) <- c('celltype',
                                  'condition')
heatmap.3(t(as.matrix(pathway_up_df_t24ht8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(11,11), main = 'upregulated pathways t24h vs t8w')
# now down
pathway_down_df_t24ht8w <- pathway_down_df[, colnames(pathway_down_df)[grep('t24ht8w', colnames(pathway_down_df))]]
pathway_down_df_t24ht8w_top_10 <- get_top_pathways(pathway_down_df_t24ht8w, 10, T)
heatmap.3(t(as.matrix(pathway_down_df_t24ht8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(t24ht8w_vary_colors_m), margins=c(11,11), main = 'downregulated pathways t24h vs t8w')

# UT vs t8w
pathway_up_df_utt8w <- pathway_up_df[, colnames(pathway_up_df)[grep('UTt8w', colnames(pathway_up_df))]]
pathway_up_df_utt8w_top_10 <- get_top_pathways(pathway_up_df_utt8w, 10, T)
utt8w_vary_colors_cond <- rep(c(cc[['UTt8w']]), times = 6)
utt8w_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
utt8w_vary_colors_m <- cbind(utt8w_vary_colors_ct, utt8w_vary_colors_cond)
colnames(utt8w_vary_colors_m) <- c('celltype',
                                     'condition')
heatmap.3(t(as.matrix(pathway_up_df_utt8w_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utt8w_vary_colors_m), margins=c(11,11), main = 'upregulated pathways HC vs t8w')

# UT vs Baseline
pathway_up_df_utbaseline <- pathway_up_df[, colnames(pathway_up_df)[grep('UTBaseline', colnames(pathway_up_df))]]
pathway_up_df_utbaseline_top_10 <- get_top_pathways(pathway_up_df_utbaseline, 10, T)
utbaseline_vary_colors_cond <- rep(c(cc[['UTBaseline']]), times = 6)
utbaseline_vary_colors_ct <- c(rep(cc[['B']], times=1),rep(cc[['CD4T']], times=1),rep(cc[['CD8T']], times=1),rep(cc[['DC']], times=1),rep(cc[['monocyte']], times=1),rep(cc[['NK']], times=1))
utbaseline_vary_colors_m <- cbind(utbaseline_vary_colors_ct, utbaseline_vary_colors_cond)
colnames(utbaseline_vary_colors_m) <- c('celltype',
                                   'condition')
heatmap.3(t(as.matrix(pathway_up_df_utbaseline_top_10)),
          col=(brewer.pal(10,"RdBu")), RowSideColors = t(utbaseline_vary_colors_m), margins=c(11,11), main = 'upregulated pathways HC vs t0')



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
pathways <- read.table('/data/scRNA/pathways/ReactomePathways.tsv', sep='\t')
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
meta.data <- read.table('/data/cardiology/metadata/cardio.integrated_meta.data.tsv', sep='\t', header=T, row.names=1)
plot_de_vs_cell_type_numbers(mast_meta_output_loc_lfc01, meta.data, plot_separately = T, proportion = T)
plot_de_vs_cell_type_numbers(mast_meta_output_loc_lfc01, meta.data, plot_separately = F, proportion = F)


# specifically for monocytes, check the number of DE genes and the fractional differences of their sub-celltype populations
plot_de_number_vs_subcell_population(mast_meta_output_loc_lfc01, 'monocyte', 'mono 1', meta.data, timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', make_absolute=F)

plot_de_number_vs_subcell_population(mast_meta_output_loc_lfc01, 'monocyte', 'mono 1', meta.data[meta.data$chem=='V2',], timepoints=c('UTBaseline', 'UTt24h', 'UTt8w', "Baselinet24h", "Baselinet8w", "t24ht8w"), cell_type_column_higherres='cell_type', cell_type_column_lowerres='cell_type_lowerres', make_absolute=F)

