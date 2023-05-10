########## Multiple testing correction for cell type composition differences for the elife method #############
####################################################################################
# Elife method obtained externally from https://elifesciences.org/articles/43882   #
####################################################################################

####################################################################################
# libraries                                                                        #
####################################################################################
library(Seurat)

####################################################################################
# Perform post-hoc Benjamin-Hochberg correction on elife p-values with W=0.1      #
####################################################################################
# Perform post-hoc Benjamin-Hochberg correction on elife p-values with W=0.1 #
allbaselinet24h <- read.table("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_baseline_t24h_20200710.tsv")
allbaselinet24h <- allbaselinet24h[c("0.1","0.05"), ]
for(colname in colnames(allbaselinet24h)){
  allbaselinet24h[[paste(colname,"corrected", sep=".")]] <- p.adjust(allbaselinet24h[[colname]], n=6, method = c("hochberg"))
}

tables <- list.files("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/", include.dirs = F)
for(table in tables){
  try({
    elife <- read.table(paste("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/", table, sep=""), header=T, sep='\t')
    elife <- elife[c("0.1","0.05"), ]
    for(colname in colnames(elife)){
      elife[[paste(colname,"hochberg", sep=".")]] <- p.adjust(elife[[colname]], n=6, method = c("hochberg"))
      elife[[paste(colname,"bonferroni", sep=".")]] <- p.adjust(elife[[colname]], n=6, method = c("bonferroni"))
      elife[[paste(colname,"holm", sep=".")]] <- p.adjust(elife[[colname]], n=6, method = c("holm"))
    }
    write.table(elife, paste("/home/umcg-ivanblokland/stemi/cell_type_composition/", table, "cor.tsv", sep=''), sep='\t')
  })
}
