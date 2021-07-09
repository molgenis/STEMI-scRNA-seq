# 6-7-2021, IV van Blokland, Roy Oelen

# Load librahaaaries
library(lme4)
library(afex)
library(lmerTest)
# ugly tidyverse library, bleh
library(tidyverse)
library(psycho)

## Loading formulas
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

# Correlating protein expression with peak CK-MB values
correlate_protein_expression_and_clinicalvars <- function(protein_expression_matrix, clinicalvar_matrix, variables_to_correlate=c("peak_ck_mb"), method='spearman', mtc_method="bonferroni"){
  # results are put into a df
  results_df <- NULL
  # get the common participants
  genes <- setdiff(colnames(protein_expression_matrix), c("id", "timepoint"))
  # check each timepoint
  for(timepoint in unique(protein_expression_matrix$timepoint)) {
    # grab the participants for the timepoint
    parts_timepoint <- unique(protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, "id"])
    # subset to specific timepoint
    protein_expression_matrix_timepoint <- protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, ]
    # check each gene
    for(gene in genes){
      # check each clinical variable
      for(clinical_variable in intersect(colnames(clinicalvar_matrix), variables_to_correlate)) {
        # get the participants that also have the variable we are interested in
        common_participants <- intersect(parts_timepoint, clinicalvar_matrix[!is.na(clinicalvar_matrix[[clinical_variable]]), "record_id"])
        # get expression data
        clinvar_values <- as.vector(unlist(clinicalvar_matrix[match(common_participants, clinicalvar_matrix$record_id), clinical_variable]))
        protein_expression <- as.vector(unlist(protein_expression_matrix_timepoint[match(common_participants, protein_expression_matrix_timepoint$id), gene]))
        # calculate the correlation
        p.value <- NA
        estimate <- NA
        correlation <- NA
        try({
          # inside a try block so if we can't correlate, we'll still do the rest
          correlation <- cor.test(clinvar_values, protein_expression, method = method)
          # get p-value of the correlation
          p.value <- correlation$p.value
          estimate <- correlation$estimate[["rho"]]
        })
        # results are put in a row
        row <- data.frame(timepoint=c(timepoint), gene=c(gene), clinical_variable=c(clinical_variable), correlation=c(estimate), p_value=c(p.value), stringsAsFactors = F)
        # add the row to the table
        if(is.null(results_df)) {
          results_df <- row
        }
        else{
          results_df <- rbind(results_df, row)
        }
      }
    }
  }
  results_df[[mtc_method]] <- p.adjust(results_df$p_value, method = mtc_method)
  return(results_df)
}

# Plotting protein expression with peak CK-MB values
plot_protein_expression_and_clinicalvars <- function(protein_expression_matrix, clinicalvar_matrix, variables_to_correlate=c("peak_ck_mb"), method='spearman', mtc_method="bonferroni"){
  # results are put into lists
  plots_per_timepoint <- list()
  # get the common participants
  genes <- setdiff(colnames(protein_expression_matrix), c("id", "timepoint"))
  # check each timepoint
  for(timepoint in unique(protein_expression_matrix$timepoint)) {
    # put the results per gene in a list
    plots_per_gene <- list()
    # grab the participants for the timepoint
    parts_timepoint <- unique(protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, "id"])
    # subset to specific timepoint
    protein_expression_matrix_timepoint <- protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, ]
    # check each gene
    for(gene in genes){
      # save per clinical variable
      plots_per_variable <- list()
      
      # check each clinical variable
      for(clinical_variable in intersect(colnames(clinicalvar_matrix), variables_to_correlate)) {
        # get the participants that also have the variable we are interested in
        common_participants <- intersect(parts_timepoint, clinicalvar_matrix[!is.na(clinicalvar_matrix[[clinical_variable]]), "record_id"])
        # get expression data
        clinvar_values <- as.vector(unlist(clinicalvar_matrix[match(common_participants, clinicalvar_matrix$record_id), clinical_variable]))
        protein_expression <- as.vector(unlist(protein_expression_matrix_timepoint[match(common_participants, protein_expression_matrix_timepoint$id), gene]))
        # put the data into a plot frame
        plot_data <- data.frame(x=protein_expression, y=clinvar_values)
        # sort
        plot_data <- plot_data[order(plot_data$x), ]
        # try to plot
        try({
          plotted <- ggplot(data=plot_data, mapping = aes(x=x, y=y)) + geom_point() + geom_smooth(method = "lm") + xlab("protein expression") + ylab(clinical_variable)
          # put in the list
          plots_per_variable[[clinical_variable]] <- plotted
        })
        # results are put in a row
      }
      plots_per_gene[[gene]] <- plots_per_variable
    }
    plots_per_timepoint[[timepoint]] <- plots_per_gene
  }
  return(plots_per_timepoint)
}

# calculating peak CK-MB values of every participant
add_peak_ck_mb <- function(clinvar_object){
  # subset to the ckmb variables
  ckmb <- clinvar[, grep('ck_mb', colnames(clinvar_object))]
  # change the NA variables to -1
  ckmb[is.na(ckmb)] <- -1
  # get the maximum
  max_ckmb <- apply(ckmb, 1, max)
  # replace a -1 with NA, this means there was no real variable
  max_ckmb[max_ckmb == -1] <- NA
  # paste onto the existing object
  clinvar_object$peak_ck_mb <- max_ckmb
  return(clinvar_object)
}

# performing mixed models to test whether peak CK-MB values im STEMI participants are associated to protein expression
# starting with the simplest model
perform_mixed_model_1 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id"){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(protein, "~", clinvar_column, "+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # perform the actual analysis
    result <- lmer(data = this_protein_data, formula = formula)
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_2 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id"){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(protein, "~", clinvar_column, "+", timepoint_column,"+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # perform the actual analysis
    result <- lmer(data = this_protein_data, formula = formula)
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_3 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id"){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(protein, "~", clinvar_column, "+", timepoint_column, "+", clinvar_column, ":", timepoint_column, "+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # perform the actual analysis
    result <- lmer(data = this_protein_data, formula = formula)
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# putting all the mixed model outcomes in a pretty table
lmm_list_to_table <- function(lmm_output){
  resulting_table <- NULL
  # check each output
  for(protein in names(lmm_output)){
    # get from the summary data
    results_pretty <- data.frame(summary(lmm_output[[protein]])$coefficients)
    # add  the rownames as a variable
    results_pretty$variable <- rownames(results_pretty)
    # set in nicer order
    results_pretty <- results_pretty[, c(ncol(results_pretty), 1:ncol(results_pretty) -1)]
    # grab the resulting colnames
    result_columns <- colnames(results_pretty)
    # add the protein as a column
    results_pretty$protein <- protein
    # set it in a nicer order
    results_pretty <- results_pretty[, c("protein", result_columns)]
    # add to the entire frame
    if(is.null(resulting_table)){
      resulting_table <- results_pretty
    }
    else {
      resulting_table <- rbind(resulting_table, results_pretty)
    }
  }
  # remove the rownames
  rownames(resulting_table) <- NULL
  return(resulting_table)
}

# permutation function
# perform every model on evey protein
model_per_protein <- function(protein_data, clinvar_data, timepoint_column="timepoint", clinvar_column="peak_ck_mb", id_column="id", clinvar_id_column="record_id"){
  resulting_list <- list()
  # for every protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # make a list for protein to test 
    list_protein <- list()
    # perform every model 
    lmm_model1 <- perform_mixed_model_1(protein_data=protein_data[, c(protein, id_column, timepoint_column)], clinvar_data=clinvar_data, clinvar_column=clinvar_column, timepoint_column=timepoint_column, id_column=id_column, clinvar_id_column=clinvar_id_column)
    lmm_model2 <- perform_mixed_model_2(protein_data=protein_data[, c(protein, id_column, timepoint_column)], clinvar_data=clinvar_data, clinvar_column=clinvar_column, timepoint_column=timepoint_column, id_column=id_column, clinvar_id_column=clinvar_id_column)
    lmm_model3 <- perform_mixed_model_3(protein_data=protein_data[, c(protein, id_column, timepoint_column)], clinvar_data=clinvar_data, clinvar_column=clinvar_column, timepoint_column=timepoint_column, id_column=id_column, clinvar_id_column=clinvar_id_column)
    # perform F-tests between models
    # interpretation: if p-value = 0.13 --> Model 1 has a greater variance and so model 2 is the best fit
    ftest_model1_model2 <- var.test(residuals(lmm_model1[[protein]]),residuals(lmm_model2[[protein]]), alternative="greater")    # if p-value of the F-test is >0.05 --> save outcomes of model 1
    ideal_model <- NA
    model_used <- NA
    # check 1 vs 2
    if(ftest_model1_model2$p.value > 0.05){
      ideal_model <- lmm_model1
      model_used <- 1
    }
    else{
      ideal_model <- lmm_model2
      model_used <- 2
    }
    # check vs. 3
    ftest_vs_model3 <- var.test(residuals(ideal_model[[protein]]),residuals(lmm_model3[[protein]]), alternative="greater")
    # overwrite if this model is better
    if(ftest_vs_model3$p.value < 0.05){
      ideal_model <- lmm_model3
      model_used <- 3
    }
    # save the result to a list
    list_protein[["model"]] <- ideal_model
    list_protein[["version"]] <- model_used
    resulting_list[[protein]] <- list_protein
  }
  return(resulting_list)
}
  
# permutation function for performing the best model fit several times
perform_permutations <- function(protein_data, clinvar_data, number_of_permutations, timepoint_column="timepoint", clinvar_column="peak_ck_mb", id_column="id", clinvar_id_column="record_id"){
  permutations_list <- list()
  # loop through all permutations
  for(i in 1:number_of_permutations){
    # 
    permuted_ids <- clinvar_data[sample(1:nrow(clinvar_data), nrow(clinvar_data), replace = F), clinvar_id_column]
    # copy the true table
    permuted_clinvar_data <- clinvar_data
    # set the permuted IDs instead
    permuted_clinvar_data[[clinvar_id_column]] <- permuted_ids
    # run function choosing model of best fit again
    permutated_model <- model_per_protein(protein_data, clinvar_data = permuted_clinvar_data)
    # add permitation to a list
    permutations_list[[i]] <- permutated_model
  }
  return(permutations_list)
}

# function to calculate FDR
perform_protein_analysis_clinvars <- function(protein_data, clinvar_data, number_of_permutations, mean_permuted_ps=T, timepoint_column="timepoint", clinvar_column="peak_ck_mb", id_column="id", clinvar_id_column="record_id"){
  # do non-permuted analysis
  real_models <- model_per_protein(protein_data, clinvar_data)
  # do permuted analysis
  permutated_models <- perform_permutations(protein_data, clinvar_data, number_of_permutations)
  # we need to store the p-values somewhere
  permuted_ps_per_protein <- list()
  # get permuted p-values per protein, check each permutation
  for(i in 1:number_of_permutations){
    # get that specific permutation
    list_this_round <- permutated_models[[i]]
    # check each protein
    for(protein in names(list_this_round)){
      # get the specific p value
      pval <- summary(list_this_round[[protein]][["model"]][[protein]])$coefficients[clinvar_column, "Pr(>|t|)"]
      # add the p value to the list
      if(protein %in% names(permuted_ps_per_protein)){
        # if this was not the first round, we can just add the p-value
        permuted_ps_per_protein[[protein]] <-c(permuted_ps_per_protein[[protein]], pval)
      }
      else{
        # if this is the first round, we need to a new vector
        permuted_ps_per_protein[[protein]] <- c(pval)
      }
    }
  }
  # add all these p values together
  permuted_ps <- c()
  for(protein in names(permuted_ps_per_protein)){
    # grab that specific protein
    ps_this_permuted_protein <- permuted_ps_per_protein[[protein]]
    # add mean to all ps if requested
    if(mean_permuted_ps){
      permuted_ps <- c(permuted_ps, mean(ps_this_permuted_protein))
    }
    # otherwise just add all
    else{
      permuted_ps <- c(permuted_ps, ps_this_permuted_protein)
    }
  }
  # store an fdr per protein
  fdr_per_protein <- list()
  # compare p-values between real and permutated
  for(protein in names(real_models)){
    # get each protein from the actual result and look at the fraction of permuted Ps, that has a better or equal P-value
    actual_p <- summary(real_models[[protein]][["model"]][[protein]])$coefficients[clinvar_column, "Pr(>|t|)"]
    # check against the permuted ones
    fdr <- length(permuted_ps[permuted_ps < actual_p])/length(permuted_ps)
    # add to the list
    fdr_per_protein[[protein]] <- fdr
  }
  return(fdr_per_protein)
}  

####################
#    Main Code     #
####################
# location of the clinical variables
clinvar.loc <- '/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/REDCap/SingleCellSequencing_20210114.tsv'
# location of the protein variables
olinkid_to_uid_loc <- '/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/20210701_royoelen/cardiology/olinkid_to_uniprotid.tsv'
uniprotid_to_gs_loc <- '/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/20210701_royoelen/cardiology/uniprot_to_genesymbol.tsv'
gs_to_ens_loc <- '/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/20210701_royoelen/cardiology/eQTL_mapping/features_v3.tsv'
olink_loc <- '/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/20210701_royoelen/cardiology/20200442_Groot_NPX-QC_format_fixed.tsv'

# load some tables
gs_to_ens <- read.table(gs_to_ens_loc, sep = '\t', header = F)
uniprotid_to_gs <- read.table(uniprotid_to_gs_loc, sep = '\t', header = T)
olinkid_to_uid <- read.table(olinkid_to_uid_loc, sep = '\t', header = T)

# read with rownames instead
olink <- read.table(olink_loc, sep = '\t', header = T, row.names = 1)
olink_plottable <- olink_to_plottable_table(olink)
# remove the test samples
olink_plottable <- olink_plottable[olink_plottable$timepoint %in% c('Baseline', 't24h', 't8w'), ]

# read the clinvar table
clinvar <- read.table(clinvar.loc, header = T, sep = '\t', dec = '.')
# calculate the peak ck_mb
clinvar <- add_peak_ck_mb(clinvar)
# add log of peak and baseline ckmb
clinvar$log_peak_ck_mb <- log(clinvar$peak_ck_mb)
clinvar$log_ck_mb <- log(clinvar$ck_mb)

# calculate correlation
corgenes <- correlate_protein_expression_and_clinicalvars(olink_plottable, clinvar, variables_to_correlate = c("peak_ck_mb"), mtc_method = "fdr")
corgenes$uniprotid <- olinkid_to_uid[match(corgenes$gene, olinkid_to_uid$OlinkID), "Uniprot.ID"]
corgenes$gs <- uniprotid_to_gs[match(corgenes$uniprotid, uniprotid_to_gs$From), "To"]
write.table(corgenes, file="/Users/irene/Documents/Geneeskunde/MD-PhD/sahacRNA-seq_data/Olink/20210701_corrclinvarstoproteins/20210701_peakckmbtoprotein.tsv", row.names= F, sep="\t")

# plotting significantly correlated proteins with peak CK-MB values
cor_ntprobnp24h <- NA
cor_ntprobnp24h <- olink_plottable[which(olink_plottable$timepoint == "t24h"),(olink_plottable=="OID00634")]
plot(clinvar$peak_ck_mb, cor_ntprobnp24h)

best_plots_eva <- plot_protein_expression_and_clinicalvars(olink_plottable, clinvar, variables_to_correlate = c("peak_ck_mb"))
best_plots_eva <- plot_protein_expression_and_clinicalvars(olink_plottable, clinvar, variables_to_correlate = c("peak_ck_mb"))
best_plots_eva[["t24h"]][["OID00634"]][["peak_ck_mb"]]

# import dataset with peak CK-MB values
olink_peakckmb <- read.table('/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/Olink_results/20210702_olink_plottable2.tsv', header = T, sep = '\t', dec = '.')

# performing LMM for model 1 and saving as model 1
lmm_protein_peakckmb_1 <- perform_mixed_model_1(olink_plottable, clinvar)
lmm_table_1 <- lmm_list_to_table(lmm_protein_peakckmb_1)
write.table(lmm_table_1, file="/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/Olink/20210701_corrclinvarstoproteins/20210705_lmm_1.tsv", row.names= F, sep="\t")

# performing LMM for model 2 and saving as model 2
lmm_protein_peakckmb_2 <- perform_mixed_model_2(olink_plottable, clinvar)
lmm_table_2 <- lmm_list_to_table(lmm_protein_peakckmb_2)
write.table(lmm_table_2, file="/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/Olink/20210701_corrclinvarstoproteins/20210705_lmm_2.tsv", row.names= F, sep="\t")

# performing LMM for model 3 and saving as model 3
lmm_protein_peakckmb_3 <- perform_mixed_model_3(olink_plottable, clinvar)
lmm_table_3 <- lmm_list_to_table(lmm_protein_peakckmb_3)
write.table(lmm_table_3, file="/Users/irene/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/Olink/20210701_corrclinvarstoproteins/20210705_lmm_3.tsv", row.names= F, sep="\t")

# F-test
# testing protein OID00616
# performing a one-sided F-test to compare the variance of the residuals of Method 1 and Method 2. We test whether Method 1 is better fit than 2, and so our H0 is: the residuals of method 1 are greater than that of method 2.
var.test(residuals(lmm_protein_peakckmb_1[["OID00616"]]),residuals(lmm_protein_peakckmb_2[["OID00616"]]), alternative="greater")
# p-value = 0.001 --> Model 1 has greater variance and so model 2 is better fitting

# performing a one-sided F-test to compare the variance of the residuals of method 2 and method 3
lm1 <- var.test(residuals(lmm_protein_peakckmb_2[["OID00616"]]),residuals(lmm_protein_peakckmb_3[["OID00616"]]), alternative="greater")
# p-value = 0.13 --> Model 2 does not have a greater variance and so model 2 is the best fit

# choosing the model of best fit
best_model <- model_per_protein(protein_data=olink_plottable, clinvar_data=clinvar)

# performing permutations
permutations_model <- perform_permutations(protein_data=olink_plottable, clinvar_data=clinvar, number_of_permutations=3) 
# perform 100 permutations so we can detect a small enough p-value
perform_protein_analysis_clinvars(protein_data=olink_plottable, clinvar_data=clinvar, number_of_permutations=100)

# significant proteins associated with peak CK-MB values 
# $OID00634; p-value=0.0326087 --> IL1RL1. Model used: 2 (protein ~ clinvar_column + timepoint_column + (1|id_column)
# $OID00600; p-value=0.04347826 --> MPO. Model used: 2 (protein ~ clinvar_column + timepoint_column + (1|id_column)
# $OID00590; p-value=0.01086957 --> TFPI. Model used:2 (protein ~ clinvar_column + timepoint_column + (1|id_column)
# $OID00131; p-value=0 --> NT-proBNP. Model used: 2 (protein ~ clinvar_column + timepoint_column + (1|id_column)
# $OID00571; p-value=0.02173913 --> TNFRSF11B. Model used:2 (protein ~ clinvar_column + timepoint_column + (1|id_column)

# plotting the significant proteins
olink_plottable$peakckmb <- clinvar[match(olink_plottable$id, clinvar$record_id), "peak_ck_mb"]
plot1 <- ggplot(data=olink_plottable, mapping=aes(x = peakckmb, y = OID00634, colour = timepoint)) + geom_point() + geom_smooth(method="lm") + theme_classic() + labs(y="IL1RL1")
plot2 <- ggplot(data=olink_plottable, mapping=aes(x = peakckmb, y = OID00600, colour = timepoint)) + geom_point() + geom_smooth(method="lm") + theme_classic() + labs(y="MPO")
plot3 <- ggplot(data=olink_plottable, mapping=aes(x = peakckmb, y = OID00590, colour = timepoint)) + geom_point() + geom_smooth(method="lm") + theme_classic() + labs(y="TFPI")
plot4 <- ggplot(data=olink_plottable, mapping=aes(x = peakckmb, y = OID00131, colour = timepoint)) + geom_point() + geom_smooth(method="lm") + theme_classic() + labs(y="NT-proBNP")
plot5 <- ggplot(data=olink_plottable, mapping=aes(x = peakckmb, y = OID00571, colour = timepoint)) + geom_point() + geom_smooth(method="lm") + theme_classic() + labs(y="TNFRSF11B")

ggsave("~/Documents/Geneeskunde/MD-PhD/scRNA-seq_data/Olink/20210615_propcorr/20210708_peakckmb_corr_protein_plots/plot_tnfrsf11b_peakckmb.pdf", dpi=600, width=20, height=20, units="cm")
plot4+plot3+plot5+plot1+plot2+plot_layout(guides = "collect")

# generating plots
# NT-proBNP
summary(best_model$OID00131$model$OID00131)
confint(best_model$OID00571$model$OID00571)


