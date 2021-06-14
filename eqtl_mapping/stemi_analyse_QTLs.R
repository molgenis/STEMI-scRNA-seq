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




# load object
#cardio.stemi <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio.stemi.20210611.combatcorrected.rds')
# genotype location
snps_loc<-'/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
# load genotype data
snps <- fread(snps_loc, header=T)
# locations of features
features_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/'
