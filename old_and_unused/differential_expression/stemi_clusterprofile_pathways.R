#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen, Paola Pibiri
# Name: stemi_clusterprofile_pathways.R
# Function: run cluster profiler on the DE output
############################################################################################################################

####################
# libraries        #
####################
library(clusterProfiler)
library(organism, character.only = TRUE)
library(enrichplot)
library(pathview)
library(ggnewscale)
library(ggplot2)


####################
# Functions        #
####################

do_cluster_profiler <- function(mast.output.loc, pathway.output.loc, cell.types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), condition.combinations=c('Baselinet24h', 'Baselinet8w'), log.fc.column='metafc',nPerm=10000, minGSSize=3, maxGSSize=800, pvalueCutoff=0.05, pAdjustMethod='BH'){
  # set the organism
  kegg.organism <- "hsa"
  organism <- 'org.Hs.eg.db'
  # backup the working directory
  prev.wd <- getwd()
  # check each cell type
  for(cell.type in cell.types){
    # check each condition combination
    for(condition.combination in condition.combinations){
      # get the full path to the DE file
      de.output.loc <- paste(mast.output.loc, cell.type, condition.combination, '.tsv', sep = '')
      # read the output
      de.output <- read.table(de.output.loc, header = T, row.names = 1, sep = '\t', stringsAsFactors = F)
      # get the entrez IDs for the gene symbols
      ids <- bitr(rownames(de.output), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)
      # get the IDs which are not duplicated, and are not NA
      ids.not.duplicated <- ids[!is.na(duplicated(ids[['ENTREZID']])) & !duplicated(ids[['SYMBOL']]) & !duplicated(ids[['ENTREZID']]), ]
      # subset the DE output to not have duplicates or NA
      de.output <- de.output[rownames(de.output) %in% ids.not.duplicated[['SYMBOL']], ]
      # add the entrez IDs as a column
      de.output[['entrezid']] <- ids[match(rownames(de.output), ids[['SYMBOL']]), 'ENTREZID']
      # order the table by the log fold change
      de.output <- de.output[order(de.output[[log.fc.column]], decreasing = T), ]
      # grab the log FCs
      log.fcs <- de.output[[log.fc.column]]
      # and set the entrez IDs as names for a named vector
      names(log.fcs) <- de.output[['entrezid']]
      # run clusterprofiler
      cc.result <- gseKEGG(geneList     = log.fcs,
                     organism     = kegg.organism,
                     nPerm        = nPerm,
                     minGSSize    = minGSSize,
                     maxGSSize    = maxGSSize,
                     pvalueCutoff = pvalueCutoff,
                     pAdjustMethod = pAdjustMethod,
                     keyType       = "ncbi-geneid")

      # extract the resulting table
      cc.table <- cc.result@result
      # check if there is a result
      if(nrow(cc.table) > 0){
        # paste the output path together
        cc.table.loc <- paste(pathway.output.loc, 'table/', cell.type, '_', condition.combination, '.tsv', sep = '')
        if(!dir.exists(paste(pathway.output.loc, 'table/', sep = ''))){
          # create it if it does not
          dir.create(paste(pathway.output.loc, 'table/', sep = ''), recursive = T)
        }
        # write the result
        write.table(cc.table, cc.table.loc, sep = '\t', quote = F, row.names = F, col.names = T)

        # create the dotplot
        cc.dotplot <- dotplot(cc.result, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
        # paste the output path together
        cc.dotplot.folder <- paste(pathway.output.loc, 'dotplot/', sep = '')
        cc.dotplot.file <- paste(cell.type, '_', condition.combination, '.pdf', sep = '')
        if(!dir.exists(cc.dotplot.folder)){
          # create it if it does not
          dir.create(cc.dotplot.folder, recursive = T)
        }
        # save the dotplot
        ggsave(filename = cc.dotplot.file, path = cc.dotplot.folder, plot = cc.dotplot, width = 10, height = 10)

        # get the similarity matrix
        cc.simmatrix <- pairwise_termsim(cc.result)
        # create the emmaplot
        cc.emmaplot <- emapplot(cc.simmatrix)
        # paste the output path together
        cc.emmaplot.folder <- paste(pathway.output.loc, 'emmaplot/', sep = '')
        cc.emmaplot.file <- paste(cell.type, '_', condition.combination, '.pdf', sep = '')
        if(!dir.exists(cc.emmaplot.folder)){
          # create it if it does not
          dir.create(cc.emmaplot.folder, recursive = T)
        }
        # save the emmaplot
        ggsave(filename = cc.emmaplot.file, path = cc.emmaplot.folder, plot = cc.emmaplot, width = 10, height = 10)

        # now we need to set the working directory for the pathview plots
        pathview.plot.folder <- paste(pathway.output.loc, 'pathview/', cell.type, '_', condition.combination, '/', sep = '')
        # check if the directory exists
        if(!dir.exists(pathview.plot.folder)){
          # create it if it does not
          dir.create(pathview.plot.folder, recursive = T)
        }
        # set the working directory to be this directory
        setwd(pathview.plot.folder)
        # check each path id
        for(path.id in unique(cc.table[['ID']])){
          # for some reason pathwview may error
          tryCatch({
            pathview(gene.data=log.fcs, pathway.id = path.id, species = kegg.organism)
          }, error=function(error_condition) {
            print(paste("pathview failed for:", cell.type, condition.combination, path.id, '->', error_condition))
          })
        }
      }
    }
  }
  # put back the working directory
  setwd(prev.wd)
}


####################
# Main Code        #
####################

# the path of the mast meta output
mast.meta.output.loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/differential_expression/MAST/stemi_meta_paired_lores_lfc01minpct01ncountrna_20210301/rna/'
# the path of the pathway output
pathway.output.loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/pathway_enrichment/kegg/MAST/stemi_meta_paired_lores_lfc01minpct01ncountrna_20210301/rna/'
# and the specific plots/tables
pathway.output.tables.loc <- paste(pathway.output.loc, 'table/', sep = '')
pathway.output.dotplots.loc <- paste(pathway.output.loc, 'dotplot/', sep = '')
pathway.output.emmaplots.loc <- paste(pathway.output.loc, 'emmaplot/', sep = '')
pathway.output.pathviewplots.loc <- paste(pathway.output.loc, 'pathview/', sep = '')
# do the profiler
do_cluster_profiler(mast.output.loc = mast.meta.output.loc, pathway.output.loc = pathway.output.loc, condition.combinations = c('UTBaseline', 'Baselinet24h', 'Baselinet8w', 't24ht8w'), pvalueCutoff = 0.05)
