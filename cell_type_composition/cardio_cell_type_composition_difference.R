####################################################################################
# libraries                                                                        #
####################################################################################
library(Seurat)
library(ggplot2)
library(parallel)
library(future.apply)
library(lme4)
library(scdney)
library(dunn.test)

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
  cond_names<-dimnames(obs.table)[[1]];
  clus_names<-dimnames(obs.table)[[2]];
  
  cond<-rep(cond_names, apply(obs.table,1,sum) );
  clus_list<-apply(obs.table,1,function(x){
    rep(clus_names, x );
  })
  clus<-base::do.call(c,clus_list);
  #table(cond,clus);
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
get_cell_counts <- function(metadata){
  # init count matrix
  cell_counts <- matrix(nrow=length(unique(metadata$timepoint.final)), ncol = length(unique(metadata$cell_type_lowerres)), dimnames = list(unique(metadata$timepoint.final), unique(metadata$cell_type_lowerres)))
  # check each condition
  for(condition in unique(metadata$timepoint.final)){
    total_cells_condition <- nrow(metadata[metadata$timepoint.final == condition, ])
    # check each cell type
    for(cell_type in unique(metadata$cell_type_lowerres)){
      cells_type_condition <- nrow(metadata[metadata$timepoint.final == condition & metadata$cell_type_lowerres == cell_type, ])
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
      res = two.class.test(obs.counts, null_dist, cond.control=cond1, cond.treatment=cond2,to.plot=F)
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
  results <- mclapply(pairs, function(pair,cell_counts, null_distributions){
    # grab the two conditions
    cond1 <- pair[1]
    cond2 <- pair[2]
    # init result
    res.table = c()
    # go through the err probs
    for(err_prob in ps){
      print(paste(cond1, cond2, err_prob))
      # obtain the null dist for the error prop that was supplied
      null_dist <- null_distributions[[as.character(err_prob)]]
      # get the result for this pair and err prob
      res = two.class.test(obs.counts, null_dist, cond.control=cond1, cond.treatment=cond2,to.plot=F)
      # add the result
      res.table <- rbind(res.table, res)
    }
    # set the names of this table as the err probs that we used
    rownames(res.table) <- as.character(ps)
    # add this table to the list by the name of the pair
    return(res.table)
  }
  ,cell_counts, null_distributions, mc.cores = ncores)
  names(results) <- names(pairs)
  return (results)
}

perform_wilcoxon_rank_sum <- function(metadata, p.adjust.method='bonferroni'){
#  # let's try wilcoxon rank sums
#  w.rank.sum.cts <- list()
#  for(cell_type in unique(metadata$cell_type_lowerres)){
#    timepointsvector <- c()
#    proportions <- c()
#    # check for timepoints
#    for(timepoint in unique(metadata$timepoint.final)){
#      # check the participants
#      for(participant in unique(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type, ]$assignment.final)){
#        # get the number of cells for the complete combination
#        ct_cond_part_number <- nrow(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type & metadata$assignment.final == participant, ])
#        # get the total number of cells regardsless of cell type
#        cond_part_number <- nrow(metadata[metadata$timepoint.final == timepoint & metadata$assignment.final == participant, ])
#        # get proportion
#        proportion <- ct_cond_part_number/cond_part_number
#        # add to proportions
#        proportions <- c(proportions, proportion) # I know this is horribly inefficient, will fix this
#      }
#      # add the number of timepoints
#      timepointsvector <- c(timepointsvector, rep(timepoint, length(unique(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type, ]$assignment.final))))
#    }
#    # do the test
#    w.rank.sum <- pairwise.wilcox.test(proportions, timepointsvector)
#    # add the results
#    w.rank.sum.cts[[cell_type]] <- w.rank.sum
#  }
  w.rank.sum.cts <- perform_categorical_test(metadata, use_proportion = T, test.use = 'wilcoxon.rank.sum', p.adjust.method = p.adjust.method)
  return(w.rank.sum.cts)
}

perform_kruskal_wallace <- function(metadata, use_proportion=T, p.adjust.method='bonferroni'){
  k.wallace.cts <- perform_categorical_test(metadata, use_proportion = use_proportion, test.use = 'kruskal.wallace', p.adjust.method = p.adjust.method)
  return(k.wallace.cts)
}

perform_dunn_test <- function(metadata, use_proportion=T, p.adjust.method='bonferroni'){
  dunn.cts <- perform_categorical_test(metadata, use_proportion = use_proportion, test.use = 'dunn', p.adjust.method = p.adjust.method)
  return(dunn.cts)
}

perform_categorical_test <- function(metadata, use_proportion=T, test.use='kruskal.wallace', p.adjust.method='bonferroni'){
  # let's try wilcoxon rank sums
  tests.cts <- list()
  for(cell_type in unique(metadata$cell_type_lowerres)){
    timepointsvector <- c()
    values <- c()
    # check for timepoints
    for(timepoint in unique(metadata$timepoint.final)){
      # check the participants
      for(participant in unique(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type, ]$assignment.final)){
        # get the number of cells for the complete combination
        ct_cond_part_number <- nrow(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type & metadata$assignment.final == participant, ])
        if(use_proportion){
          # get the total number of cells regardsless of cell type
          cond_part_number <- nrow(metadata[metadata$timepoint.final == timepoint & metadata$assignment.final == participant, ])
          # get proportion
          value <- ct_cond_part_number/cond_part_number
          # add to proportions
          values <- c(values, value) # I know this is horribly inefficient, will fix this
        }
        else{
          values <- c(values, ct_cond_part_number) # I know this is horribly inefficient, will fix this
        }
        
      }
      # add the number of timepoints
      timepointsvector <- c(timepointsvector, rep(timepoint, length(unique(metadata[metadata$timepoint.final == timepoint & metadata$cell_type_lowerres == cell_type, ]$assignment.final))))
    }
    # do the test
    test <- NULL
    if(test.use == 'dunn'){
      test <- dunn.test(values, timepointsvector, kw = F, method = c(p.adjust.method))
    }
    else if(test.use == 'kruskal.wallace'){
      test <- k.wallace <- kruskal.test(values, timepointsvector, p.adjust.method = p.adjust.method)
    }
    else if(test.use == 'wilcoxon.rank.sum'){
      test <- pairwise.wilcox.test(proportions, timepointsvector, p.adjust.method = p.adjust.method)
    }
    # add the results
    tests.cts[[cell_type]] <- test
  }
  return(tests.cts)
}

perform_fisher_exact <- function(cell_counts){
  fisher.result.cts <- list()
  # let's try the fisher exact test as well
  for(cell_type in colnames(cell_counts)){
    # get the number of cells for the cell type per condition
    nr_of_cell_type <- cell_counts[,cell_type]
    # get the number of cells per condition
    total_number_of_cells <- rowSums(cell_counts)
    # get the non-cell type cells per condition
    nr_of_cells_not_cell_type <- total_number_of_cells - nr_of_cell_type
    # create frame
    contingency_table <- data.frame(nr_of_cell_type, nr_of_cells_not_cell_type, row.names = rownames(cell_counts))
    colnames(contingency_table) <- c(cell_type, 'other')
    # doing fisher.exact on all conditions at once is apparently to heavy, so we need to do it pairwise
    result_table <- data.frame(matrix(NA, nrow = length(rownames(contingency_table)), ncol = length(rownames(contingency_table))))
    colnames(result_table) <- rownames(cell_counts)
    rownames(result_table) <- rownames(cell_counts)
    # do the comparisons
    for(i in 1:length(rownames(contingency_table))){
      i_start <- i+1
      if(i_start <= length(rownames(contingency_table))){
        for(i2 in i_start:length(rownames(contingency_table))){
          cond1 <- rownames(contingency_table)[i]
          cond2 <- rownames(contingency_table)[i2]
          # subset the contingency table
          contingency_subtable <- contingency_table[c(cond1, cond2),]
          # perform the test
          fisher.result <- fisher.test(t(contingency_subtable))
          # set in the result table
          result_table[cond1, cond2] <- fisher.result$p.value
          result_table[cond2, cond1] <- fisher.result$p.value
          }
      }
    }
    # add result
    fisher.result.cts[[cell_type]] <- result_table
  }
  return(fisher.result.cts)
}

fitGLM <- function(res, condition, subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE){
  
  fit_random <- list()
  fit_fixed <- list()
  for(i in 1:ncol(res$nstar)){
    # idx <- indexes_list[[i]]
    if(verbose){
      if(i%%10==0){
        print(paste("fitting GLM...", i))
      }
    }
    
    
    glm_df <-  cbind(res$info[,1:2], res$nstar[,i])
    
    # glm_df <- melt(glm_df)
    colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
    glm_df$cond <- condition
    
    if(subject_effect){
      if(pairwise){
        
        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df, family = poisson(link=log))
        if(!fixed_only){
          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond +  cellTypes:cond + (1 | subject ),
                                         data = glm_df, family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }else{
        fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df, family = poisson(link=log))
        if(!fixed_only){
          
          fit_random[[i]] <- lme4::glmer(cell_count ~ cellTypes + cond + cellTypes:cond + (1 | subject ),
                                         data = glm_df, family = poisson(link=log),
                                         control = glmerControl(nAGQ = 0L))
        }
      }
    }else{
      fit_fixed[[i]] <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond, data = glm_df,
                                   family = poisson(link=log))
    }
    
  }
  
  if(!subject_effect){
    fixed_only = TRUE
  }
  
  if(!fixed_only){
    pool_res_random = mice::pool(fit_random)
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_random = pool_res_random,
                pool_res_fixed = pool_res_fixed,
                fit_random = fit_random,
                fit_fixed = fit_fixed))
  }else{
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_fixed = pool_res_fixed,
                fit_fixed = fit_fixed))
  }
  
}

fitGLM_mt <- function(res, condition, subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE){
  fit_fixed <- future_apply(res$nstar, 2, function(x){  
    print('a thread has started work')
    glm_df <-  cbind(res$info[,1:2], x)
    colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
    glm_df$cond <- condition
    if(subject_effect){
      if(pairwise){
        print('a thread has started a fixed glm')
        fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df,
                                     family = poisson(link=log))
      }else{
        print('a thread has started a fixed glm')
        fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                     data = glm_df, family = poisson(link=log))
      }
    }else{
      print('a thread has started a fixed glm')
      fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond, data = glm_df, family = poisson(link=log))
    }
    print('a thread has finished work')
    return(fit_fixed_i)
  })
  fit_random <- list()
  if(!subject_effect){
    fixed_only = TRUE
  }
  else{
    if(pairwise){
      fit_random <- future_apply(res$nstar, 2, function(x){
        print('a thread has started a paired glm')
        glm_df <-  cbind(res$info[,1:2], x)
        colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
        glm_df$cond <- condition
        fit_random_i <- lme4::glmer(cell_count ~ cellTypes + cond + cellTypes:cond + (1 | subject ), data = glm_df, family = poisson(link=log), control = glmerControl(nAGQ = 0L))
        print('a thread has finished work')
        return(fit_random_i)
      })
    }
    else{
      fit_random <- future_apply(res$nstar, 2, function(x){
        print('a thread has started a paired glm')
        glm_df <-  cbind(res$info[,1:2], x)
        colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
        glm_df$cond <- condition
        fit_random_i <- lme4::glmer(cell_count ~ cellTypes + cond + cellTypes:cond + (1 | subject ), data = glm_df, family = poisson(link=log), control = glmerControl(nAGQ = 0L))
        print('a thread has finished work')
        return(fit_random_i)
      })
    }
  }
  
  if(!fixed_only){
    pool_res_random = mice::pool(fit_random)
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_random = pool_res_random,
                pool_res_fixed = pool_res_fixed,
                fit_random = fit_random,
                fit_fixed = fit_fixed))
  }else{
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_fixed = pool_res_fixed,
                fit_fixed = fit_fixed))
  }
  
}


fitGLM_mc <- function(res, condition, subject_effect = TRUE, pairwise = TRUE, fixed_only = FALSE, verbose = TRUE, mc.cores = 1){
  res_list <- list()
  for(i in 1:ncol(res$nstar)){
    res_list[[i]] <- res$nstar[,i]
  }
  fit_fixed <- mclapply(res_list, function(x, info){
    print('a thread has started work')
    glm_df <-  cbind(info, x)
    colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
    glm_df$cond <- condition
    if(subject_effect){
      if(pairwise){
        print('a thread has started a fixed glm')
        fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                  data = glm_df,
                                  family = poisson(link=log))
      }else{
        print('a thread has started a fixed glm')
        fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond + subject,
                                  data = glm_df, family = poisson(link=log))
      }
    }else{
      print('a thread has started a fixed glm')
      fit_fixed_i <- stats::glm(cell_count ~ cellTypes + cond +  cellTypes:cond, data = glm_df, family = poisson(link=log))
    }
    print('a thread has finished work')
    return(fit_fixed_i)
  }, res$info[,1:2], mc.cores = mc.cores, mc.preschedule = F)
  fit_random <- list()
  if(!subject_effect){
    fixed_only = TRUE
  }
  else{
    if(pairwise){
      fit_random <- mclapply(res_list, function(x, info){
        print('a thread has started a paired glm')
        glm_df <-  cbind(info, x)
        colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
        glm_df$cond <- condition
        fit_random_i <- lme4::glmer(cell_count ~ cellTypes + cond + cellTypes:cond + (1 | subject ), data = glm_df, family = poisson(link=log), control = glmerControl(nAGQ = 0L))
        print('a thread has finished work')
        return(fit_random_i)
      }, res$info[,1:2], mc.cores = mc.cores, mc.preschedule = F)
    }
    else{
      fit_random <- mclapply(res_list, function(x, info){
        print('a thread has started a paired glm')
        glm_df <-  cbind(info, x)
        colnames(glm_df) <- c("cellTypes", "subject", "cell_count")
        glm_df$cond <- condition
        fit_random_i <- lme4::glmer(cell_count ~ cellTypes + cond + cellTypes:cond + (1 | subject ), data = glm_df, family = poisson(link=log), control = glmerControl(nAGQ = 0L))
        print('a thread has finished work')
        return(fit_random_i)
      }, res$info[,1:2], mc.cores = mc.cores, mc.preschedule = F)
    }
  }
  
  if(!fixed_only){
    pool_res_random = mice::pool(fit_random)
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_random = pool_res_random,
                pool_res_fixed = pool_res_fixed,
                fit_random = fit_random,
                fit_fixed = fit_fixed))
  }else{
    pool_res_fixed = mice::pool(fit_fixed)
    return(list(pool_res_fixed = pool_res_fixed,
                fit_fixed = fit_fixed))
  }
  
}


####################################################################################
# main code                                                                        #
####################################################################################

# read Seurat file
cardio.integrated <- readRDS('/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/objects/cardio_integrated_20pcs_ctTmonoTc17Th17_20200622.rds')
# grab metadata
metadata <- cardio.integrated@meta.data
# grab all the cell counts
cell_count_all <- get_cell_counts(metadata)
# get the null distributions
null_dist_all <- get_null_distributions_mt(cell_count_all)
# create the pairs we which to test
ut_baseline <- c('UT', 'Baseline')
baseline_t8w <- c('Baseline', 't8w')
ut_t8w <- c('UT', 't8w')
pairs <- list()
pairs[['ut_baseline']] <- ut_baseline
pairs[['baseline_t8w']] <- baseline_t8w
pairs[['ut_t8w']] <- ut_t8w
# get the results
diff_all <- test_two_class_mt(cell_count_all, null_dist_all, pairs)
# write results
write.table(diff_all[['ut_baseline']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_ut_baseline_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_all[['ut_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_ut_t8w_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_all[['baseline_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/all_baseline_t8w_20200705.tsv', sep = '\t', row.names=T)

metadata_v2 <- metadata[metadata$chem == 'V2', ]
cell_count_v2 <- get_cell_counts(metadata_v2)
# get the null distributions
null_dist_v2 <- get_null_distributions_mt(cell_count_v2)
# get the results
diff_v2 <- test_two_class_mt(cell_count_v2, null_dist_v2, pairs)
# write results
write.table(diff_v2[['ut_baseline']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v2_ut_baseline_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['ut_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v2_ut_t8w_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_v2[['baseline_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v2_baseline_t8w_20200705.tsv', sep = '\t', row.names=T)

metadata_v3 <- metadata[metadata$chem == 'V3', ]
cell_count_v3 <- get_cell_counts(metadata_v3)
# get the null distributions
null_dist_v3 <- get_null_distributions_mt(cell_count_v3)
# get the results
diff_v3 <- test_two_class_mt(cell_count_v3, null_dist_v3, pairs)
# write results
write.table(diff_v3[['ut_baseline']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v3_ut_baseline_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['ut_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v3_ut_t8w_20200705.tsv', sep = '\t', row.names=T)
write.table(diff_v3[['baseline_t8w']], '/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/elife_2020/v3_baseline_t8w_20200705.tsv', sep = '\t', row.names=T)




# let's try the chi2 methods as wel
chisq_all <- chisq.test(cell_count_all)
chisq_v2 <- chisq.test(cell_count_v2)
chisq_v3 <- chisq.test(cell_count_v3)

# fisher exact as well
fisher_all <- perform_fisher_exact(metadata)
fisher_v2 <- perform_fisher_exact(metadata_v2)
fisher_v3 <- perform_fisher_exact(metadata_v3)

# kruskal-wallace too
k.wallace_all <- perform_kruskal_wallace(metadata)
k.wallace_v2 <- perform_kruskal_wallace(metadata_v2)
k.wallace_v3 <- perform_kruskal_wallace(metadata_v3)

# do a dunn post-hoc test
dunn_all <- perform_dunn_test(metadata)
dunn_v2 <- perform_dunn_test(metadata_v2)
dunn_v3 <- perform_dunn_test(metadata_v3)


# let's try the wilcoxon rank sums test
w.rank.sum.cts.all <- perform_wilcoxon_rank_sum(metadata)
w.rank.sum.cts.v2 <- perform_wilcoxon_rank_sum(metadata_v2)
w.rank.sum.cts.v3 <- perform_wilcoxon_rank_sum(metadata_v3)
# check the cell types
for(ct in names(w.rank.sum.cts.all)){
  write.table(w.rank.sum.cts.all[[ct]]$p.value, paste('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/wilcoxon_ranked_sum/all_', ct, '.tsv', sep = ''), sep = '\t', row.names = T)
  write.table(w.rank.sum.cts.v2[[ct]]$p.value, paste('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/wilcoxon_ranked_sum/v2_', ct, '.tsv', sep = ''), sep = '\t', row.names = T)
  write.table(w.rank.sum.cts.v3[[ct]]$p.value, paste('/groups/umcg-wijmenga/scr01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/wilcoxon_ranked_sum/v3_', ct, '.tsv', sep = ''), sep = '\t', row.names = T)
}


# now trying the scDC method
# grab metadata
metadata <- cardio.integrated@meta.data
#exprsMat <- sim$sim_exprsMat
exprsMat <- cardio.integrated@assays$SCT@data
#subject <- sim$sim_subject
subject <- paste(cardio.integrated@meta.data$assignment.final, cardio.integrated@meta.data$timepoint.final, sep = '.')
#cellTypes <- sim$sim_cellTypes
cellTypes <- cardio.integrated@meta.data$cell_type_lowerres
#cond <- sim$sim_cond
cond <- cardio.integrated@meta.data$timepoint.final
# try to run scDC
res_scDC_noClust <- scDC_noClustering(cellTypes, subject, calCI = TRUE, 
                                      calCI_method = c("percentile", "BCa", "multinom"), ncores = 12, nboot = 10000)
# we need to make a condition vector
conds <- c()
for(participant in sort(unique(metadata$assignment.final))){
  # get the conditions for this participant
  for(condition in sort(unique(metadata[metadata$assignment.final == participant, ]$timepoint.final))){
    conds <- c(conds, rep(condition, length(unique(metadata$cell_type_lowerres))))
  }
}

# view proportions
barplotCI(res_scDC_noClust, conds)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/all_fixed_20200705_barplot.png', dpi = 600, height = 10, width = 10)
# view density
densityCI(res_scDC_noClust, conds)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/all_fixed_20200705_density.png', dpi = 600, height = 10, width = 10)
# try fitting a GLM
#plan(multiprocess)
res_GLM <- fitGLM_mc(res_scDC_noClust, conds, pairwise = F, mc.cores = 12)
# save the results
write.table(summary(res_GLM$pool_res_fixed), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/all_fixed_20200705.tsv', row.names = F)
write.table(summary(res_GLM$pool_res_random), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/all_random_20200705.tsv', row.names = F)

# now for v2
cardio.integrated.v2 <- subset(cardio.integrated, subset = chem == 'V2')
# grab metadata.v2
metadata.v2 <- cardio.integrated.v2@meta.data
#metadata.v2 <- sim$sim_metadata.v2
exprsMat.v2 <- cardio.integrated.v2@assays$SCT@data
#metadata.v2 <- sim$sim_metadata.v2
subject.v2 <- paste(cardio.integrated.v2@meta.data$assignment.final, cardio.integrated.v2@meta.data$timepoint.final, sep = '.')
#cellTypes.v2 <- sim$sim_cellTypes.v2
cellTypes.v2 <- cardio.integrated.v2@meta.data$cell_type_lowerres
#cond.v2 <- sim$sim_cond.v2
cond.v2 <- cardio.integrated.v2@meta.data$timepoint.final
# try to run scDC
res_scDC_noClust.v2 <- scDC_noClustering(cellTypes.v2, subject.v2, calCI = TRUE, 
                                         calCI_method = c("percentile", "BCa", "multinom"), ncores = 12, nboot = 10000)
# we need to make a condition vector
conds.v2s <- c()
for(participant in sort(unique(metadata.v2$assignment.final))){
  # get the conditions for this participant
  for(condition in sort(unique(metadata.v2[metadata.v2$assignment.final == participant, ]$timepoint.final))){
    conds.v2s <- c(conds.v2s, rep(condition, length(unique(metadata.v2$cell_type_lowerres))))
  }
}

# view proportions
barplotCI(res_scDC_noClust.v2, conds.v2s)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v2_fixed_20200705_barplot.png', dpi = 600, height = 10, width = 10)
# view density
densityCI(res_scDC_noClust.v2, conds.v2s)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v2_fixed_20200705_density.png', dpi = 600, height = 10, width = 10)
# try fitting a GLM
res_GLM.v2 <- fitGLM_mc(res_scDC_noClust.v2, conds.v2s, pairwise = F, mc.cores = 12)
# save the results
write.table(summary(res_GLM.v2$pool_res_fixed), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v2_fixed_20200705.tsv', row.names = F)
write.table(summary(res_GLM.v2$pool_res_random), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v2_random_20200705.tsv', row.names = F)

# now for v3
cardio.integrated.v3 <- subset(cardio.integrated, subset = chem == 'V3')
# grab metadata.v3
metadata.v3 <- cardio.integrated.v3@meta.data
#metadata.v3 <- sim$sim_metadata.v3
exprsMat.v3 <- cardio.integrated.v3@assays$SCT@data
#metadata.v3 <- sim$sim_metadata.v3
subject.v3 <- paste(cardio.integrated.v3@meta.data$assignment.final, cardio.integrated.v3@meta.data$timepoint.final, sep = '.')
#cellTypes.v3 <- sim$sim_cellTypes.v3
cellTypes.v3 <- cardio.integrated.v3@meta.data$cell_type_lowerres
#cond.v3 <- sim$sim_cond.v3
cond.v3 <- cardio.integrated.v3@meta.data$timepoint.final
# try to run scDC
res_scDC_noClust.v3 <- scDC_noClustering(cellTypes.v3, subject.v3, calCI = TRUE, 
                                         calCI_method = c("percentile", "BCa", "multinom"), ncores = 12, nboot = 10000)
# we need to make a condition vector
conds.v3s <- c()
for(participant in sort(unique(metadata.v3$assignment.final))){
  # get the conditions for this participant
  for(condition in sort(unique(metadata.v3[metadata.v3$assignment.final == participant, ]$timepoint.final))){
    conds.v3s <- c(conds.v3s, rep(condition, length(unique(metadata.v3$cell_type_lowerres))))
  }
}

# view proportions
barplotCI(res_scDC_noClust.v3, conds.v3s)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v3_fixed_20200705_barplot.png', dpi = 600, height = 10, width = 10)
# view density
densityCI(res_scDC_noClust.v3, conds.v3s)
ggsave('/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v3_fixed_20200705_density.png', dpi = 600, height = 10, width = 10)
# try fitting a GLM
res_GLM.v3 <- fitGLM_mc(res_scDC_noClust.v3, conds.v3s, pairwise = F, mc.cores = 12)
# save the results
write.table(summary(res_GLM.v3$pool_res_fixed), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v3_fixed_20200705.tsv', row.names = F)
write.table(summary(res_GLM.v3$pool_res_random), '/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cell_type_composition/scdney/v3_random_20200705.tsv', row.names = F)
