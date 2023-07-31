#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Marc-Jan Bonder, Roy Oelen
# Name: stemi_limix_mtc.R
# Function: do multiple testing correction on the limix output
############################################################################################################################

####################
# libraries        #
####################

library(optparse) # NOT IN CONTAINER
library(qvalue)

####################
# Functions        #
####################

####################
# Main Code        #
####################

# make command line options
option_list <- list(
  make_option(c("-d", "--directory"), type="character", default='./', 
              help="input directory", metavar="character"),
  make_option(c("-c", "--cutoff"), type="numeric", default=0.1,
              help="FDR significance cutoff", metavar="numeric")
)

# initialize optparser
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# set the directory
directory <- opt[['directory']]
fdrCutoff <- opt[['cutoff']]

# read the full output file
qtls <- read.table(paste(directory, "/qtl_results_all.txt", sep = ''), header = T, sep = '\t', stringsAsFactors = F)
# create a column of the combination of the variant and the feature
qtls["QTL"] = paste(qtls$snp_id,qtls$feature_id,sep=":")
# order by the emperical, then variant-feature p values
qtls = qtls[order(qtls$empirical_feature_p_value,qtls$p_value),]
# for variants that are in multiple chunks, take the first entry, which should be the most significant one
qtls = qtls[which(!duplicated(qtls$QTL)),]

# remove variants that behave strange. (strict)
qtls = qtls[which(qtls$alpha_param>=0.25),]
qtls = qtls[which(qtls$alpha_param<=4.0),]
# only use p values which are not zero or one
qtls = qtls[intersect(which(qtls$empirical_feature_p_value>=0),which(qtls$empirical_feature_p_value<=1)),]

# take the first QTL for each feature, which should be the top one, as we sorted them before
topQtls = qtls[which(!duplicated(qtls$feature_id)),]
# calculate a global FDR using qvaules
topQtls["Global_FDR"] <- qvalue::qvalue(topQtls$empirical_feature_p_value)$qvalue

# all effects with the global FDR
write.table(topQtls,paste0(directory, "/top_qtl_results_all_fdr.txt"),quote=F, sep="\t",col.names = T, row.names=F)
# siignificant top effects
write.table(topQtls[which(topQtls$Global_FDR<=fdrCutoff),],paste0(directory, "/top_qtl_results_all_filtered_fdr",fdrCutoff,".txt"),quote=F, sep="\t",col.names = T, row.names=F)

# get the max emperical p where the qvalue FDR is below the FDR cutoff
sigCutoff = max(topQtls$empirical_feature_p_value[which(topQtls$Global_FDR<=fdrCutoff)])
# subset to QTLs an emperical p value below this cutoff
sigQtls <- qtls[which(qtls$empirical_feature_p_value<=sigCutoff),]
# write significant subset
write.table(sigQtls ,paste0(directory, "/qtl_results_all_filtered_fdr",fdrCutoff,".txt"),quote=F, sep="\t",col.names = T, row.names=F)