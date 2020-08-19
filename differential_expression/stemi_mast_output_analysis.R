####################
# libraries        #
####################

library(RColorBrewer)
library(MetaVolcanoR)
library(ggplot2)
library(data.table)


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

get_top_vary_genes <- function(de_table, use_tp=T, use_pathogen=T, use_ct=T, sd_cutoff=0.5, use_dynamic_sd=F, top_so_many=10, must_be_positive_once=F, pathogens=c("CA", "MTB", "PA"), timepoints=c("3h", "24h"), cell_types=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  top_vary_de <- c()
  cols_to_loop <- NULL
  # grab the appriate grep
  if(use_tp & use_pathogen){
    # we want a combo of pathogen and timepoints, so 3hCA for example
    cols_to_loop <- paste(rep(timepoints, each = length(pathogens)), pathogens, sep = "")
  }
  else if(use_pathogen & use_ct){
    # cell type and pathogen, so monocyte3hCA and monocyte24hCA for example
    cols_to_loop <- paste(rep(cell_types, each = length(pathogens)), pathogens, sep = ".*")
  }
  else if(use_tp & use_ct){
    # cell type at a timepoint, so monocyte3hCA and monocyte3hPA and monocyte3hMTB for example
    cols_to_loop <- paste(rep(cell_types, each = length(timepoints)), timepoints, sep = ".*")
  }
  else if(use_pathogen){
    cols_to_loop <- pathogens
  }
  else if(use_tp){
    cols_to_loop <- timepoints
  }
  else if(use_ct){
    cols_to_loop <- cell_types
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
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# get the location of the pathways
v2_pathway_up_output_loc <- '/data/cardiology/pathways/sigs_pos/v2_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
v2_pathway_up_df <- get_pathway_table(v2_pathway_up_output_loc, append = '_sigpos.txt')
v2_pathway_up_df[v2_pathway_up_df==0] <- 350
# get the df limited by top pathways of upregulated genes
v2_pathway_up_df_top_3 <- get_top_pathways(v2_pathway_up_df, 3, T)

heatmap.3(t(as.matrix(v2_pathway_up_df_top_3)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# get the location of the pathways
pathway_down_output_loc <- '/data/cardiology/pathways/sigs_neg/meta_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
pathway_down_df <- get_pathway_table(pathway_down_output_loc, append = '_pathwayneg.txt')
pathway_down_df[pathway_down_df==0] <- 400
# get the df limited by top pathways of upregulated genes
pathway_down_df_top_3 <- get_top_pathways(pathway_down_df, 3, T)

heatmap.3(t(as.matrix(pathway_down_df_top_3)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

# get the location of the pathways
v2_pathway_down_output_loc <- '/data/cardiology/pathways/sigs_neg/v2_paired_lores_lfc01minpct01_20200707_ensid/rna/'
# write the combined pathway file
v2_pathway_down_df <- get_pathway_table(v2_pathway_down_output_loc, append = '_signeg.txt')
v2_pathway_down_df[v2_pathway_down_df==0] <- 400
# get the df limited by top pathways of upregulated genes
v2_pathway_down_df_top_3 <- get_top_pathways(v2_pathway_down_df, 3, T)

heatmap.3(t(as.matrix(v2_pathway_down_df_top_3)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_m), margins=c(15,10))

