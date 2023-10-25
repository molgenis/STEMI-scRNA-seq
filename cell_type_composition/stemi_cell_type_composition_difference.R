#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_cell_type_composition_difference.R
# Function: perform gene set enrichment with REACTOME and enrichr
############################################################################################################################


####################################################################################
# libraries                                                                        #
####################################################################################
library(Seurat)
library(ggplot2)
library(parallel)
library(future.apply)

####################################################################################
# functions                                                                        #
####################################################################################

####################################################################################
# v <obtained externally from https://elifesciences.org/articles/43882 > v #
####################################################################################

# Differential proportion analysis of single-cell populations

###
# Simulation of a data set with different number of cells and their cell types
###
simCells<-function(n=10^7, prop=c(0.6, 0.3, 0.05, 0.03, 0.01, 0.005, 0.005)){
  ct<-paste("Cell_Type_",1:length(prop),sep="");
  nct<-round(n*prop);
  lab<-rep(ct,nct);	
  return(lab);
}

###
# Simulate a data-set with some provided error rate to simulate experimental/biological noise
###
simCellsNoisy<-function(n=10^7, prop, error = 0.05){
  # Create a 'noisy' distribution for sampling cells
  prop.noisy = unlist(lapply(prop, function(x) 
  {ifelse(sample(c(1:2), size=1) == 1, x+x*error, x-x*error)}))
  prop.noisy = prop.noisy/sum(prop.noisy)
  
  # Sample as before but with a modified distribution
  ct<-paste("Cell_Type_",1:length(prop.noisy),sep="");
  nct<-round(n*prop.noisy);
  lab<-rep(ct,nct);	
  return(lab);
}

###
# Subsampling of cells
###
subCells<-function(cells, n=round(length(cells)/10)){
  sample(cells,n);
}



###
# Create the condition and cluster variables from the observed contingency table
###
makeVariables<-function(obs.table){
  # grab the condition names, 
  cond_names<-dimnames(obs.table)[[1]];
  # [1] "t24h"     "Baseline" "t8w"      "UT"
  
  # grab the cluster/celltype names
  clus_names<-dimnames(obs.table)[[2]];
  # [1] "CD4T"             "Treg"             "CD8T"             "monocyte"        
  # [5] "DC"               "NK"               "B"                "megakaryocyte"   
  # [9] "plasma B"         "Th17"             "hemapoietic stem" "Tc17"
  
  # for the sum of the counts of one condition, repeat the timepoints/conditions
  cond<-rep(cond_names, apply(obs.table,1,sum) );
  # chr [1:141270] "t24h" "t24h" "t24h" "t24h" "t24h" "t24h" "t24h" "t24h"
  
  # for the rows, repeat the cluster/ct names for the times the cell/clus appeared in the condition
  clus_list<-apply(obs.table,1,function(x){
    rep(clus_names, x );
  })
  # List of 4
  # $ t24h    : chr [1:36218] "CD4T" "CD4T" "CD4T" "CD4T" ...
  # $ Baseline: chr [1:33070] "CD4T" "CD4T" "CD4T" "CD4T" ...
  # $ t8w     : chr [1:32231] "CD4T" "CD4T" "CD4T" "CD4T" ...
  # $ UT      : chr [1:39751] "CD4T" "CD4T" "CD4T" "CD4T" ...
  
  # assign each cell type/cell to it's timepoint as well
  clus<-base::do.call(c,clus_list);
  # Named chr [1:141270] "CD4T" "CD4T" "CD4T" "CD4T" "CD4T" "CD4T" "CD4T" ...
  # - attr(*, "names")= chr [1:141270] "t24h1" "t24h2" "t24h3" "t24h4" ...
  
  #table(cond,clus);
  
  # return a list where we have this list of conditions and timepoints per cell
  return(list(cond=cond, clus=clus));
}

###
# Generate null distributions based on a permutation procedure
###
generateNull<-function(obs.table, n=10000, p=0.2){
  obs<-makeVariables(obs.table);
  all.exp<-array(NA,dim=c(dim(obs.table)[1],dim(obs.table)[2],n))
  dimnames(all.exp)[[1]]<-dimnames(obs.table)[[1]];
  dimnames(all.exp)[[2]]<-dimnames(obs.table)[[2]];
  dimnames(all.exp)[[3]]<-1:n;
  clus_names<-dimnames(obs.table)[[2]];
  cond_names<-dimnames(obs.table)[[1]];
  # logi [1:4, 1:12, 1:2] NA NA NA NA NA NA ...
  # - attr(*, "dimnames")=List of 3
  # ..$ : chr [1:4] "t24h" "Baseline" "t8w" "UT"
  # ..$ : chr [1:12] "CD4T" "Treg" "CD8T" "monocyte" ...
  # ..$ : chr [1:2] "1" "2"
  
  v<-makeVariables(obs.table);
  
  ### Perform permutation n times
  for(i in 1:n){
    pv<-v;
    
    if (i %% 1000 == 0){
      print(i);}
    
    ### Permutation: Randomly select p*100% of points to be re-sampled from the background 
    randInd<-sample(1:length(pv$clus),round(length(pv$clus)*p));
    pv$clus[randInd]<-sample(v$clus,length(randInd),replace=F);  
    
    this.exp<-table(pv$cond,pv$clus);
    exp.table<-obs.table;
    exp.table[1:dim(exp.table)[1],1:dim(exp.table)[2]]<-NA
    exp.table[dimnames(this.exp)[[1]],dimnames(this.exp)[[2]]]<-this.exp;
    exp.table
    all.exp[,,i]<-exp.table;
  }
  return(all.exp);
}

get_permutated_counts <- function(obs.table, p=0.2){
  obs<-makeVariables(obs.table);
  clus_names<-dimnames(obs.table)[[2]];
  cond_names<-dimnames(obs.table)[[1]];
  
  v<-makeVariables(obs.table);
  pv<-v;
  
  randInd<-sample(1:length(pv$clus),round(length(pv$clus)*p));
  pv$clus[randInd]<-sample(v$clus,length(randInd),replace=F);  
  
  this.exp<-table(pv$cond,pv$clus);
  exp.table<-obs.table;
  exp.table[1:dim(exp.table)[1],1:dim(exp.table)[2]]<-NA
  exp.table[dimnames(this.exp)[[1]],dimnames(this.exp)[[2]]]<-this.exp;
  
  return(exp.table)
}

###
# Perform sum by ignoring NA
###
sumNoNA<-function(x){
  sum(x[which(!is.na(x))])
}

###
# Perform a two-class test of significance
###
two.class.test<-function(obs.table, all.exp, cond.control="C", cond.treatment="PA",to.plot=T){
  clus_names<-dimnames(obs.table)[[2]]
  pp<-array(-1,length(clus_names));
  names(pp)<-clus_names;
  
  if(to.plot){
    par(mfrow=c(3,ceiling(length(clus_names)/3)));
  }
  for(this_clus in clus_names){
    obs.diff<-obs.table[cond.treatment,this_clus]/sum(obs.table[cond.treatment,]) - 
      obs.table[cond.control,this_clus]/sum(obs.table[cond.control,]);
    all.diff<-all.exp[cond.treatment,this_clus,]/apply(all.exp[cond.treatment,,],2,sumNoNA) -
      all.exp[cond.control,this_clus,]/apply(all.exp[cond.control,,],2,sumNoNA);
    if(to.plot){
      hist(all.diff,breaks=50,col="grey",border="grey",main=this_clus)
      abline(v=obs.diff,col="red")
    }
    p.r<-length(which(obs.diff>all.diff))/length(which(!is.na(all.diff)));
    p.l<-length(which(obs.diff<all.diff))/length(which(!is.na(all.diff)));
    pp[this_clus]<-min(p.r,p.l);
  }
  return(pp);
}

########################################################################################
# ^  </ end of obtained externally from https://elifesciences.org/articles/43882 >  ^  #
########################################################################################

# get the cell counts from the metadata object
get_cell_counts <- function(metadata, cell_type_column='cell_type_lowerres'){
  # init count matrix
  cell_counts <- matrix(nrow=length(unique(metadata$timepoint.final)), ncol = length(unique(metadata[[cell_type_column]])), dimnames = list(unique(metadata$timepoint.final), unique(metadata[[cell_type_column]])))
  # check each condition
  for(condition in unique(metadata$timepoint.final)){
    total_cells_condition <- nrow(metadata[metadata$timepoint.final == condition, ])
    # check each cell type
    for(cell_type in unique(metadata[[cell_type_column]])){
      cells_type_condition <- nrow(metadata[metadata$timepoint.final == condition & metadata[[cell_type_column]] == cell_type, ])
      # cant devide by zero and don't want to devide zero
      cells_fraction <- 0
      if(total_cells_condition > 0 & cells_type_condition > 0 ){
        cells_fraction <- cells_type_condition/total_cells_condition
      }
      # set the value
      cell_counts[condition,cell_type] <- cells_fraction
      cell_counts[condition,cell_type] <- cells_type_condition
    }
  }
  return(cell_counts)
}


# build the null distributions of the cell counts, with different number of iterations and error prob
get_null_distributions <- function(cell_counts, n=100000, ps=c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)){
  # put the null distributions in a list
  null_distributions = list()
  # work for every p to make a null prob for
  for(err_prob in ps){
    # generate null distribution for this err prob
    tip.exp <- generateNull(cell_counts, n, p=err_prob)
    null_distributions[[as.character(err_prob)]] = tip.exp
  }
  return(null_distributions)
}

# build the null distributions of the cell counts, with different number of iterations and error prob using multiple cores
get_null_distributions_mt <- function(cell_counts, n=100000, ps=c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05), mc.cores=NULL){
  # we need characters as names for the list
  ps_character <- as.character(ps)
  # check if the number of cores was supplied
  ncores <- mc.cores
  if(is.null(ncores)){
    # by default use the number of distributions we are making
    ncores <- length(ps)
  }
  # put the null distributions in a list
  null_distributions <- mclapply(ps, function(err_prob, n, cell_counts){
    # generate null distribution for this err prob
    tip.exp <- generateNull(cell_counts, n, p=err_prob)
    return(tip.exp)
  }
  ,n,cell_counts, mc.cores = ncores)
  names(null_distributions) <- ps_character
  return(null_distributions)
}

# test the cell type composition for pairs and different error prob null distributions
test_two_class <- function(cell_counts, null_distributions, pairs, ps=c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05)){
  # store the result of the pairs in a list
  results <- list()
  # check each pair
  for(pair in names(pairs)){
    # grab the two conditions
    cond1 <- pairs[[pair]][1]
    cond2 <- pairs[[pair]][2]
    # init result
    res.table = c()
    # go through the err probs
    for(err_prob in ps){
      # obtain the null dist for the error prop that was supplied
      null_dist <- null_distributions[[as.character(err_prob)]]
      # get the result for this pair and err prob
      res = two.class.test(cell_counts, null_dist, cond.control=cond1, cond.treatment=cond2,to.plot=F)
      # add the result
      res.table <- rbind(res.table, res)
    }
    # set the names of this table as the err probs that we used
    rownames(res.table) <- as.character(ps)
    # add this table to the list by the name of the pair
    results[[pair]] <- res.table
  }
  return(results)
}

# test the cell type composition for pairs and different error prob null distributions using multiple cores
test_two_class_mt <- function(cell_counts, null_distributions, pairs, ps=c(0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05), mc.cores=NULL){
  # we need characters as names for the list
  ps_character <- as.character(ps)
  # check if the number of cores was supplied
  ncores <- mc.cores
  if(is.null(ncores)){
    # by default use the number of distributions we are making
    ncores <- length(ps)
  }
  results <- mclapply(pairs, function(pair,cell_counts, null_distributions, ws){
    # grab the two conditions
    cond1 <- pair[1]
    cond2 <- pair[2]
    # init result
    res.table = c()
    # go through the err probs
    for(err_prob in ws){
      print(paste(cond1, cond2, err_prob))
      # obtain the null dist for the error prop that was supplied
      null_dist <- null_distributions[[as.character(err_prob)]]
      # get the result for this pair and err prob
      res = two.class.test(cell_counts, null_dist, cond.control=cond1, cond.treatment=cond2,to.plot=F)
      # add the result
      res.table <- rbind(res.table, res)
    }
    # set the names of this table as the err probs that we used
    rownames(res.table) <- as.character(ps)
    # add this table to the list by the name of the pair
    return(res.table)
  }
  ,cell_counts, null_distributions, ps, mc.cores = ncores)
  names(results) <- names(pairs)
  return (results)
}


####################################################################################
# main code                                                                        #
####################################################################################

# read Seurat file
cardio.integrated <- readRDS('/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/seurat/cardio.integrated.20210301.rds')
# grab metadata
metadata <- cardio.integrated@meta.data
# grab all the cell counts
cell_count_all <- get_cell_counts(metadata)
# get the null distributions
# null_dist_all <- get_null_distributions_mt(cell_count_all)
# create the pairs we which to test
ut_baseline <- c('UT', 'Baseline')
ut_t24h <- c('UT', 't24h')
ut_t8w <- c('UT', 't8w')
baseline_t24h <- c('Baseline', 't24h')
baseline_t8w <- c('Baseline', 't8w')
t24h_t8w <- c('t24h', 't8w')
pairs <- list()
pairs[['ut_baseline']] <- ut_baseline
pairs[['ut_t24h']] <- ut_t24h
pairs[['ut_t8w']] <- ut_t8w
pairs[['baseline_t24h']] <- baseline_t24h
pairs[['baseline_t8w']] <- baseline_t8w
pairs[['t24h_t8w']] <- t24h_t8w

# # get the results
# diff_all <- test_two_class_mt(cell_count_all, null_dist_all, pairs)
# # write results
# write.table(diff_all[['ut_baseline']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_ut_baseline_20210301.tsv', sep = '\t', row.names=T)
# write.table(diff_all[['ut_t24h']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_ut_t24h_20210301.tsv', sep = '\t', row.names=T)
# write.table(diff_all[['ut_t8w']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_ut_t8w_20210301.tsv', sep = '\t', row.names=T)
# write.table(diff_all[['baseline_t24h']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_baseline_t24h_20210301.tsv', sep = '\t', row.names=T)
# write.table(diff_all[['baseline_t8w']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_baseline_t8w_20210301.tsv', sep = '\t', row.names=T)
# write.table(diff_all[['t24h_t8w']], '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_t24h_t8w_20210301.tsv', sep = '\t', row.names=T)


metadata_v2 <- metadata[metadata$chem == 'V2', ]
cell_count_v2 <- get_cell_counts(metadata_v2)
# get the null distributions
null_dist_v2 <- get_null_distributions_mt(cell_count_v2, ps=c(0.05))
# get the results
diff_v2 <- test_two_class_mt(cell_count_v2, null_dist_v2, pairs, ps=c(0.05))
# write results
write.table(diff_v2[['ut_baseline']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_ut_baseline_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['ut_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_ut_t8w_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t8w_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['baseline_t24h']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t24h_20210301.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t8w_20210301.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['t24h_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_t24h_t8w_20210301.tsv', sep = '\t', row.names=T)


metadata_v3 <- metadata[metadata$chem == 'V3', ]
cell_count_v3 <- get_cell_counts(metadata_v3)
# get the null distributions
null_dist_v3 <- get_null_distributions_mt(cell_count_v3, ps=c(0.05))
# get the results
diff_v3 <- test_two_class_mt(cell_count_v3, null_dist_v3, pairs, ps=c(0.05))
# write results
write.table(diff_v3[['ut_baseline']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_ut_baseline_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['ut_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_ut_t8w_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t8w_20201209.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['baseline_t24h']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t24h_20210301.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t8w_20210301.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['t24h_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_t24h_t8w_20210301.tsv', sep = '\t', row.names=T)

# for more granular cell types
cell_count_hires_v2 <- get_cell_counts(metadata_v2, cell_type_column = 'cell_type')
# get the null distributions
null_dist_hires_v2 <- get_null_distributions_mt(cell_count_hires_v2, ps=c(0.05))
# get the results
diff_hires_v2 <- test_two_class_mt(cell_count_hires_v2, null_dist_hires_v2, pairs, ps=c(0.05))
# write results
write.table(diff_hires_v2[['ut_baseline']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_ut_baseline_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v2[['ut_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_ut_t8w_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v2[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t8w_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v2[['baseline_t24h']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t24h_20210301_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v2[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_baseline_t8w_20210301_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v2[['t24h_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v2_t24h_t8w_20210301_hires.tsv', sep = '\t', row.names=T)

# for more granular cell types
cell_count_hires_v3 <- get_cell_counts(metadata_v3, cell_type_column = 'cell_type')
# get the null distributions
null_dist_hires_v3 <- get_null_distributions_mt(cell_count_hires_v3, ps=c(0.05))
# get the results
diff_hires_v3 <- test_two_class_mt(cell_count_hires_v3, null_dist_hires_v3, pairs, ps=c(0.05))
# write results
write.table(diff_hires_v3[['ut_baseline']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_ut_baseline_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v3[['ut_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_ut_t8w_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v3[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t8w_20201209_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v3[['baseline_t24h']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t24h_20210301_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v3[['baseline_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_baseline_t8w_20210301_hires.tsv', sep = '\t', row.names=T)
write.table(diff_hires_v3[['t24h_t8w']], '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_type_composition/elife_2020/v3_t24h_t8w_20210301_hires.tsv', sep = '\t', row.names=T)
