#
# stemi_protein_to_clinvars.R
#
# try to explain peak_ck_mb, which is a proxy for heart damage, by protein expression, measured using the CVD3 olink panel
#
# 6-7-2021, IV van Blokland, Roy Oelen
#

#################
# libraries     #
#################
# Load librahaaaries
library(lme4)
library(lmerTest)

#################
# functions     #
#################

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
correlate_protein_expression_and_clinicalvars <- function(protein_expression_matrix, clinicalvar_matrix, variables_to_correlate=c("peak_ck_mb"), method='spearman', mtc_method="bonferroni", id_column_clinvar='record_id', id_column_protein='id', timepoint_column_protein='timepoint'){
  # results are put into a df
  results_df <- NULL
  # get the common participants
  genes <- setdiff(colnames(protein_expression_matrix), c(id_column_protein, timepoint_column_protein))
  # check each timepoint
  for(timepoint in unique(protein_expression_matrix$timepoint)) {
    # grab the participants for the timepoint
    parts_timepoint <- unique(protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, id_column_protein])
    # subset to specific timepoint
    protein_expression_matrix_timepoint <- protein_expression_matrix[protein_expression_matrix$timepoint == timepoint, ]
    # check each gene
    for(gene in genes){
      # check each clinical variable
      for(clinical_variable in intersect(colnames(clinicalvar_matrix), variables_to_correlate)) {
        # get the participants that also have the variable we are interested in
        common_participants <- intersect(parts_timepoint, clinicalvar_matrix[!is.na(clinicalvar_matrix[[clinical_variable]]), id_column_clinvar])
        # get expression data
        clinvar_values <- as.vector(unlist(clinicalvar_matrix[match(common_participants, clinicalvar_matrix[[id_column_clinvar]]), clinical_variable]))
        protein_expression <- as.vector(unlist(protein_expression_matrix_timepoint[match(common_participants, protein_expression_matrix_timepoint[[id_column_protein]]), gene]))
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
plot_protein_expression_and_clinicalvars <- function(protein_expression_matrix, clinicalvar_matrix, variables_to_correlate=c("peak_ck_mb"), method='spearman', mtc_method="bonferroni", protein_id_column='id', clinvar_id_column='record_id'){
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
perform_mixed_model_1 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F){
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
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_2 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(protein, "~", timepoint_column, "+", clinvar_column,"+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_3 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F){
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
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# Trying an alternative model
perform_mixed_model_1a <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id"){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(protein, "~", clinvar_column, "*", timepoint_column, "+", "(", timepoint_column, "|", id_column, ")", sep = "")
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


test_lmm_explain_protein <- function(protein_data, clinvar_data, protein_oid, to_fit_column, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id"){
  # subset to data for that protein
  this_protein_data <- protein_data[, c(id_column, timepoint_column, protein_oid)]
  # add the clinical variable
  this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
  # build the formula
  text_formula <- paste(protein_oid, "~", clinvar_column, "+", "(1|", id_column, ")", sep = "")
  # show what we are doing
  print(paste("doing", text_formula))
  # turn into actual formula
  formula <- as.formula(text_formula)
  # subset data to cases we have all data
  this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                           !is.na(this_protein_data[[timepoint_column]]) & 
                                           !is.na(this_protein_data[[id_column]]) & 
                                           !is.na(this_protein_data[[protein_oid]]),  ]
  # perform the actual analysis
  result <- lme4::lmer(data = this_protein_data, formula = formula)
  return(result)
}


perform_mixed_model_pred_clin_1 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~", protein, "+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_pred_clin_2 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~", timepoint_column, "+", protein,"+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_pred_clin_3 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~", protein, "+", timepoint_column, "+", clinvar_column, ":", timepoint_column, "+", "(1|", id_column, ")", sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}


perform_mixed_model_pred_clin_4 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~",  "(1|", id_column, ")", '+', protein, sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}


# making the model more complex
perform_mixed_model_pred_clin_5 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~", timepoint_column, "+", "(1|", id_column, ")", '+', protein, sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}

# making the model more complex
perform_mixed_model_pred_clin_6 <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  # store results in list for now
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # build the formula
    text_formula <- paste(clinvar_column, "~", "(1|", id_column, ")", '+', timepoint_column, '+', protein, sep = "")
    # show what we are doing
    print(paste("doing", text_formula))
    # turn into actual formula
    formula <- as.formula(text_formula)
    # subset data to cases we have all data
    this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                             !is.na(this_protein_data[[timepoint_column]]) & 
                                             !is.na(this_protein_data[[id_column]]) & 
                                             !is.na(this_protein_data[[protein]]),  ]
    # scale if requested
    if(scale){
      this_protein_data[[clinvar_column]] <- scale(this_protein_data[[clinvar_column]])
    }
    # initialize variable
    result <- NA
    if(use_lme4){
      result <- lme4::lmer(data = this_protein_data, formula = formula)
    }
    else{
      # perform the actual analysis
      result <- lmer(data = this_protein_data, formula = formula)
    }
    # save the result to a list
    result_list[[protein]] <- result
  }
  return(result_list)
}


test_lmm_explain_clinical <- function(protein_data, clinvar_data, protein_oid, to_fit_column, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", include_tp=F){
  # subset to data for that protein
  this_protein_data <- protein_data[, c(id_column, timepoint_column, protein_oid)]
  # add the clinical variable
  this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
  # build the formula
  text_formula <- paste(clinvar_column, "~", protein_oid, "+", "(1|", id_column, ")", sep = "")
  if(include_tp){
    paste(clinvar_column, "~", protein_oid, "*", timepoint_column, "+", "(", timepoint_column, "|", id_column, ")", sep = "")
  }
  # show what we are doing
  print(paste("doing", text_formula))
  # turn into actual formula
  formula <- as.formula(text_formula)
  # subset data to cases we have all data
  this_protein_data <- this_protein_data[!is.na(this_protein_data[[clinvar_column]]) & 
                                           !is.na(this_protein_data[[timepoint_column]]) & 
                                           !is.na(this_protein_data[[id_column]]) & 
                                           !is.na(this_protein_data[[protein_oid]]),  ]
  # perform the actual analysis
  result <- lme4::lmer(data = this_protein_data, formula = formula)
  return(result)
}


perform_mixed_model_pred_clin_ordered <- function(protein_data, clinvar_data, clinvar_column="peak_ck_mb", timepoint_column="timepoint", id_column="id", clinvar_id_column="record_id", use_lme4=F, scale=F){
  result_list <- list()
  # check each protein
  for(protein in setdiff(colnames(protein_data), c(timepoint_column, id_column))){
    # subset to data for that protein
    this_protein_data <- protein_data[, c(id_column, timepoint_column, protein)]
    # add the clinical variable
    this_protein_data[[clinvar_column]] <- clinvar_data[match(this_protein_data[[id_column]], clinvar[[clinvar_id_column]]), clinvar_column]
    # predict on the fixed effect
    text_formula_id <- paste(clinvar_column, "~", id_column, sep = "")
    # predict on the timepoint
    text_formula_timepoint <- paste(clinvar_column, "~", timepoint_column, sep = "")
    # predict on the protein
    text_formula_protein <- paste(clinvar_column, "~", protein, sep = "")
    # initialize the result
    result_formula_id <- NA
    result_formula_timepoint <- NA
    result_formula_protein <- NA
    if(use_lme4){
      print(paste('formula:', text_formula_id))
      result_formula_id <- lm(data = this_protein_data, formula = text_formula_id)
      print(paste('formula:', text_formula_timepoint))
      result_formula_timepoint <- lm(data = this_protein_data, formula = text_formula_timepoint)
      print(paste('formula:', text_formula_protein))
      result_formula_protein <- lm(data = this_protein_data, formula = text_formula_protein)
    }
    else{
      # perform the actual analysis
      print(paste('formula:', text_formula_id))
      result_formula_id <- lmer(data = this_protein_data, formula = text_formula_id)
      print(paste('formula:', text_formula_timepoint))
      result_formula_timepoint <- lmer(data = this_protein_data, formula = text_formula_timepoint)
      print(paste('formula:', text_formula_protein))
      result_formula_protein <- lmer(data = this_protein_data, formula = text_formula_protein)
    }
    # list of these models
    models <- list('id' = result_formula_id, 'timepoint' = result_formula_timepoint, 'protein' = result_formula_protein)
    # save the result to a list
    result_list[[protein]] <- models
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


perform_simple_linear_model <- function(protein_data, clinvar_data, clinvar_columns=c('peak_ck_mb', 'age', 'sex'), mtc_methods=c('BH', 'bonferroni'), timepoint_column='timepoint', timepoints=c('Baseline', 't24h', 't8w'), id_column="id", clinvar_id_column="record_id", clinvar_prepend='TEST_', protein_start_column_nr=3){
  # put the results in a table
  result_table <- NULL
  # for convenience sake, we will turn the timepoint columns into strings
  protein_data[[timepoint_column]] <- as.character(protein_data[[timepoint_column]])
  # if we are using a prepend
  clinvar_data[[clinvar_id_column]] <- paste(clinvar_prepend, clinvar_data[[clinvar_id_column]], sep = '')
  # we will check each timepoint
  for(timepoint in intersect(timepoints, unique(protein_data[[timepoint_column]]))){
    # subset the protein data to this timepoint
    protein_data_timepoint <- protein_data[protein_data[[timepoint_column]] == timepoint, ]
    # paste each clinical variable onto it
    protein_data_timepoint[, clinvar_columns] <- clinvar_data[match(protein_data_timepoint[[id_column]], clinvar_data[[clinvar_id_column]]), clinvar_columns]
    # check each protein, the last columns are the clinical variables we pasted onto the dataframe, so we'll skip those columns, they are not protein data
    for(protein in colnames(protein_data_timepoint)[protein_start_column_nr:(ncol(protein_data_timepoint) - length(clinvar_columns))]){
      # grab the data of the clinical variables and the protein
      protein_data_timepoint_protein <- protein_data_timepoint[, c(protein, clinvar_columns)]
      # we'll rename the protein to be 'protein', this makes grabbing the results a bit simpler
      colnames(protein_data_timepoint_protein)[[1]] <- 'protein'
      # we want only the complete entries, and get these by summing the NAs in each row, if there are no NAs, the sum is zero
      protein_data_timepoint_protein <- protein_data_timepoint_protein[apply(protein_data_timepoint_protein, 1, FUN = function(x){
        return(sum(is.na(x)))
      }) == 0, ]
      # we'll build the formula now
      formula_explanatory_variables <- paste(clinvar_columns, collapse = '+')
      formula_string <- paste('protein', '~', formula_explanatory_variables)
      formula <- formula(formula_string)
      # these values we will extract
      variables_per_clinical_variable <- c('p', 'estimate', 'se', 't')
      # make combinations of these variables for each clinical variable
      clinvar_variable_result_columns_df <- expand.grid(clinvar_columns, variables_per_clinical_variable) # creates dataframe of each combination as two columns
      clinvar_variable_result_columns <- paste(clinvar_variable_result_columns_df[[1]], clinvar_variable_result_columns_df[[2]], sep = '.') # turn the rows into the combinations we want
      # let's get all thos columns together then
      result_columns <- c('protein', 'timepoint', 'n', 'r2', 'adj_r2', clinvar_variable_result_columns)
      # create a empty result first
      result_row <- data.frame(matrix(NA, nrow=1, ncol=length(result_columns), dimnames = list(NULL, result_columns)))
      # fill in what we already know
      result_row[1, 'protein'] <- protein
      result_row[1, 'timepoint'] <- timepoint
      result_row[1, 'n'] <- nrow(protein_data_timepoint_protein)
      # we will try to do linear regression now. This might fail, but we then want to give an NA result, not exit
      try({
        linear_regression <- summary(lm(formula = formula, data = protein_data_timepoint_protein))
        # extract the shared variables
        r2 = as.numeric(linear_regression$r.squared)
        adj_r2 = as.numeric(linear_regression$adj.r.squared)
        # add to the result
        result_row[1, 'r2'] <- r2
        result_row[1, 'adj_r2'] <- adj_r2
        # check of each clinical variable the results
        for(clinvar_column in clinvar_columns){
          # get all these variables
          p = as.numeric(linear_regression$coefficients[[clinvar_column, 'Pr(>|t|)']])
          estimate = as.numeric(linear_regression$coefficients[[clinvar_column, 'Estimate']])
          se = as.numeric(linear_regression$coefficients[[clinvar_column, 'Std. Error']])
          t = as.numeric(linear_regression$coefficients[[clinvar_column, 't value']])
          # set these as variables
          result_row[1, paste(clinvar_column, 'p', sep='.')] <- p
          result_row[1, paste(clinvar_column, 'estimate', sep ='.')] <- estimate
          result_row[1, paste(clinvar_column, 'se', sep = '.')] <- se
          result_row[1, paste(clinvar_column, 't', sep = '.')] <- t
        }
      })
      # let's add this result
      if(is.null(result_table)){
        result_table <- result_row
      }
      else{
        result_table <- rbind(result_table, result_row)
      }
    }
  }
  # do MTC if requested
  for(method in mtc_methods){
    # for each p we have
    for(clinvar in clinvar_columns){
      # paste the original column name together
      orig_column <- paste(clinvar, 'p', sep = '.')
      # and the new column
      new_column <- paste(clinvar, method, sep = '.')
      # do the actual calculation
      result_table[[new_column]] <- p.adjust(result_table[[orig_column]], method = method)
    }
  }
  return(result_table)
}


calc_age_in_days <- function(date_of_birth_string, split_character='/', day_of_month_col=1, month_of_year_col=2, year_col=3){
  # split the string
  split_string <- strsplit(date_of_birth_string, split = split_character)[[1]]
  day = as.numeric(split_string[day_of_month_col])
  month = as.numeric(split_string[month_of_year_col])
  year = as.numeric(split_string[year_col])
  # create new string as input for the Date library
  standard_data_string <- paste(day, month, year, sep = '/')
  # turn into date object
  birth <- as.Date(standard_data_string, format='%d/%m/%Y')
  # get now
  now_string <- Sys.Date()
  # convert to Date as well
  now <- as.Date(now_string, format='%Y-%m-%d')
  # get the time between now and the birth
  age_days <- difftime(now, birth, units = 'days')
  return(age_days)
}


perform_simple_linear_model_protein <- function(protein_data, clinvar_data, clinvar_column='peak_ck_mb', clinvar_covars=c('age', 'sex'), mtc_methods=c('BH', 'bonferroni'), timepoint_column='timepoint', timepoints=c('Baseline', 't24h', 't8w'), id_column="id", clinvar_id_column="record_id", clinvar_prepend='TEST_', protein_start_column_nr=3){
  # put the results in a table
  result_table <- NULL
  # for convenience sake, we will turn the timepoint columns into strings
  protein_data[[timepoint_column]] <- as.character(protein_data[[timepoint_column]])
  # if we are using a prepend
  clinvar_data[[clinvar_id_column]] <- paste(clinvar_prepend, clinvar_data[[clinvar_id_column]], sep = '')
  # combine the interested variable with the covariates
  clinvar_columns <- c(clinvar_column, clinvar_covars)
  # we will check each timepoint
  for(timepoint in intersect(timepoints, unique(protein_data[[timepoint_column]]))){
    # subset the protein data to this timepoint
    protein_data_timepoint <- protein_data[protein_data[[timepoint_column]] == timepoint, ]
    # paste each clinical variable onto it
    protein_data_timepoint[, clinvar_columns] <- clinvar_data[match(protein_data_timepoint[[id_column]], clinvar_data[[clinvar_id_column]]), clinvar_columns]
    # check each protein, the last columns are the clinical variables we pasted onto the dataframe, so we'll skip those columns, they are not protein data
    for(protein in colnames(protein_data_timepoint)[protein_start_column_nr:(ncol(protein_data_timepoint) - length(clinvar_columns))]){
      # grab the data of the clinical variables and the protein
      protein_data_timepoint_protein <- protein_data_timepoint[, c(protein, clinvar_columns)]
      # we'll rename the protein to be 'protein', this makes grabbing the results a bit simpler
      colnames(protein_data_timepoint_protein)[[1]] <- 'protein'
      # we want only the complete entries, and get these by summing the NAs in each row, if there are no NAs, the sum is zero
      protein_data_timepoint_protein <- protein_data_timepoint_protein[apply(protein_data_timepoint_protein, 1, FUN = function(x){
        return(sum(is.na(x)))
      }) == 0, ]
      # we'll build the formula now
      formula_string <- paste(clinvar_column, '~', 'protein')
      # if we have covariates, we will add these as well
      if(length(clinvar_covars) > 0 ){
        formula_explanatory_variables <- paste(clinvar_covars, collapse = '+')
        formula_string <- paste(formula_string, '+', formula_explanatory_variables, sep = '')
      }
      formula <- formula(formula_string)
      # these values we will extract
      variables_per_clinical_variable <- c('p', 'estimate', 'se', 't')
      # make combinations of these variables for each covariate, and the protein
      clinvar_variable_result_columns_df <- expand.grid(c('protein', clinvar_covars), variables_per_clinical_variable) # creates dataframe of each combination as two columns
      clinvar_variable_result_columns <- paste(clinvar_variable_result_columns_df[[1]], clinvar_variable_result_columns_df[[2]], sep = '.') # turn the rows into the combinations we want
      # let's get all thos columns together then
      result_columns <- c('protein', 'timepoint', 'n', 'r2', 'adj_r2', clinvar_variable_result_columns)
      # create a empty result first
      result_row <- data.frame(matrix(NA, nrow=1, ncol=length(result_columns), dimnames = list(NULL, result_columns)))
      # fill in what we already know
      result_row[1, 'protein'] <- protein
      result_row[1, 'timepoint'] <- timepoint
      result_row[1, 'n'] <- nrow(protein_data_timepoint_protein)
      # we will try to do linear regression now. This might fail, but we then want to give an NA result, not exit
      try({
        linear_regression <- summary(lm(formula = formula, data = protein_data_timepoint_protein))
        # extract the shared variables
        r2 = as.numeric(linear_regression$r.squared)
        adj_r2 = as.numeric(linear_regression$adj.r.squared)
        # add to the result
        result_row[1, 'r2'] <- r2
        result_row[1, 'adj_r2'] <- adj_r2
        # check of each covariate, and the protein
        for(variate in c('protein', clinvar_covars)){
          # get all these variables
          p = as.numeric(linear_regression$coefficients[[variate, 'Pr(>|t|)']])
          estimate = as.numeric(linear_regression$coefficients[[variate, 'Estimate']])
          se = as.numeric(linear_regression$coefficients[[variate, 'Std. Error']])
          t = as.numeric(linear_regression$coefficients[[variate, 't value']])
          # set these as variables
          result_row[1, paste(variate, 'p', sep='.')] <- p
          result_row[1, paste(variate, 'estimate', sep ='.')] <- estimate
          result_row[1, paste(variate, 'se', sep = '.')] <- se
          result_row[1, paste(variate, 't', sep = '.')] <- t
        }
      })
      # let's add this result
      if(is.null(result_table)){
        result_table <- result_row
      }
      else{
        result_table <- rbind(result_table, result_row)
      }
    }
  }
  # do MTC if requested
  for(method in mtc_methods){
    # for each p we have
    for(variate in c('protein', clinvar_covars)){
      # paste the original column name together
      orig_column <- paste(variate, 'p', sep = '.')
      # and the new column
      new_column <- paste(variate, method, sep = '.')
      # do the actual calculation
      result_table[[new_column]] <- p.adjust(result_table[[orig_column]], method = method)
    }
  }
  # add the dependent variable for clearness sake
  result_table$dependant <- clinvar_column
  return(result_table)
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

# location of the clinical variables
clinvar.loc <- '/data/cardiology/clinical_parameters/scRNAseq_clinical_data_20200910_final_notextfields.tsv'
clinvar.loc <- '/data/cardiology/clinical_parameters/SingleCellSequencing_20210114.tsv'
# location of the protein variables
olinkid_to_uid_loc <- '/data/cardiology/olinkid_to_uniprotid.tsv'
uniprotid_to_gs_loc <- '/data/cardiology/uniprot_to_genesymbol.tsv'
gs_to_ens_loc <- '/data/cardiology/eQTL_mapping/features_v3.tsv'
olink_loc <- '/data/cardiology/20200442_Groot_NPX-QC_format_fixed.tsv'


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
clinvar <- read.table(clinvar.loc, header = T, sep = '\t', dec = '.', stringsAsFactors = F)
# calculate the peak ck_mb
clinvar <- add_peak_ck_mb(clinvar)
# add log of peak and baseline ckmb
clinvar$log_peak_ck_mb <- log(clinvar$peak_ck_mb)
clinvar$log_ck_mb <- log(clinvar$ck_mb)

# add the age
clinvar$age <- apply(clinvar, 1, function(x){
  age <- calc_age_in_days(as.character(x[['date_birth']]), split_character =  '-', year_col = 1, month_of_year_col = 2, day_of_month_col = 3)
  return(age)
})

# try to explain protein by peak_ck_mb
protein_by_peakckmb <- perform_simple_linear_model(protein_data = olink_plottable[olink_plottable[['timepoint']] == 't8w', ], clinvar_data = clinvar, clinvar_prepend='', clinvar_columns = c('log_peak_ck_mb'))
# write result
result_reg_loc <- paste('/data/cardiology/clinical_parameters', '/', 'protein_by_peakckmb.tsv', sep = '')
write.table(protein_by_peakckmb, result_reg_loc, col.names = T, row.names = F, sep = '\t')



# We first copy the olink data and put it in a metadata file
metadata <- olink_plottable[, c('id', 'timepoint')]
# Then we only take the protein data from the metadata
protein_data_log <- olink_plottable[, setdiff(colnames(olink_plottable), c('id', 'timepoint'))]
# And we log transform the proteins
protein_data_log <- log(protein_data_log)
# Add the log transformed proteins to the protein data together with the metadata
protein_data_log <- cbind(metadata, protein_data_log)
# try to explain protein by peak_ck_mb
log_protein_by_peakckmb <- perform_simple_linear_model(protein_data = protein_data_log[protein_data_log[['timepoint']] == 't8w', ], clinvar_data = clinvar, clinvar_prepend='', clinvar_columns = c('peak_ck_mb'))
result_log_reg_loc <- paste('/data/cardiology/clinical_parameters', '/', 'log_protein_by_peakckmb.tsv', sep = '')
write.table(log_protein_by_peakckmb, result_log_reg_loc, col.names = T, row.names = F, sep = '\t')


# the protein data is what you give as the variable you are interested
# olink_plottable for normal protein data
# protein_data_log for log transformed protein data

# clinvar_column is the column name of what you are trying to explain
# 'peak_ck_mb' for explaining peak_ck_mb
# 'log_peak_ck_mb' for explaining log transformed peak_ck_mb

# clinvar_covars conatins the variables you want to 'correct' for (this just means they also explain peak_ck_mb for example)

# explain peak_ck_mb by protein+age+gender
peak_ck_mb_by_protein <- perform_simple_linear_model_protein(protein_data = olink_plottable[olink_plottable[['timepoint']] == 'Baseline', ], clinvar_data = clinvar, clinvar_column='peak_ck_mb', clinvar_covars=c('age', 'gender'), clinvar_prepend='')
# explain log peak_ck_mb by protein+log_ck_mb+gender+age
log_peak_ck_mb_by_protein <- perform_simple_linear_model_protein(protein_data = olink_plottable[olink_plottable[['timepoint']] == 'Baseline', ], clinvar_data = clinvar, clinvar_column='log_peak_ck_mb', clinvar_covars=c('age', 'gender'), clinvar_prepend='')
# explain peak_ck_mb by log_protein+age+gender
peak_ck_mb_by_log_protein <- perform_simple_linear_model_protein(protein_data = protein_data_log[protein_data_log[['timepoint']] == 'Baseline', ], clinvar_data = clinvar, clinvar_column='peak_ck_mb', clinvar_covars=c('age', 'gender'), clinvar_prepend='')
# write the results to files
peak_ck_mb_by_protein_loc <- paste('/data/cardiology/clinical_parameters', '/', 'peak_ck_mb_by_protein.tsv', sep = '')
write.table(peak_ck_mb_by_protein, peak_ck_mb_by_protein_loc, sep = '\t', col.names = T, row.names = F)
log_peak_ck_mb_by_protein_loc <- paste('/data/cardiology/clinical_parameters', '/', 'log_peak_ck_mb_by_protein.tsv', sep = '')
write.table(log_peak_ck_mb_by_protein, log_peak_ck_mb_by_protein_loc, sep = '\t', col.names = T, row.names = F)
peak_ck_mb_by_log_protein_loc <- paste('/data/cardiology/clinical_parameters', '/', 'peak_ck_mb_by_log_protein.tsv', sep = '')
write.table(peak_ck_mb_by_log_protein, peak_ck_mb_by_log_protein_loc, sep = '\t', col.names = T, row.names = F)

