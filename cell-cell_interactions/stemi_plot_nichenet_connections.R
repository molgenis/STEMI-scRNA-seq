

interactions_to_circle <- function(interactions_per_ct_list){
  # get the total number of connections per cell type
  interaction_numbers <- get_interaction_numbers(interactions_per_ct_list)
  # slim it down into a three column input
  connections_slim <- slim_df_down(interaction_numbers, new_col_names = c('sender', 'receiver', 'connections'))
  # add the receiver and connections together
  connections_directionless <- combine_columns(connections_slim, c('sender', 'receiver'), 'connections')
  
  print(connections_slim)
  print(connections_directionless)
  
  # start making the plot, set height of track (cell types)
  circos.par('track.height' = 0.1)
  circos.initialize(connections_slim$sectors, x = df$x)
}

combine_columns <- function(df, cols_to_combine, variable_to_combine, na_to_zero=T){
  # init new table
  combined_df <- NULL
  if(na_to_zero){
    df[is.na(df)] <- 0
  }
  # get the possible variables
  possible_variables <- unique(c(df[[cols_to_combine[1]]], df[[cols_to_combine[2]]]))
  # check each variable
  for(possible_variable_1 in possible_variables){
    # against each other value
    for(possible_variable_2 in setdiff(possible_variables, possible_variable_1)){
      # get the combined value
      summed_value <- sum(
        as.vector(
          unlist(
            df[(df[[ cols_to_combine[1] ]] == possible_variable_1 & df[[ cols_to_combine[2] ]] == possible_variable_2) |
               (df[[ cols_to_combine[1] ]] == possible_variable_2 & df[[ cols_to_combine[2] ]] == possible_variable_1)
               , variable_to_combine]
          )
        )
      )
      # create a row
      combined_row <- data.frame(a = c(possible_variable_1), b = c(possible_variable_2), c=c(summed_value))
      # add to dataframe if requested
      if(is.null(combined_df)){
        combined_df <- combined_row
      }
      else{
        combined_df <- rbind(combined_df, combined_row)
      }
    }
  }
  return(combined_df)
}


slim_df_down <- function(df_rownames_and_colnames, new_col_names=c('x', 'y', 'z')){
  slim_df <- NULL
  # check each row
  for(row_index in 1:nrow(df_rownames_and_colnames)){
    # get the rowname, which is the first value in our input
    col1_value <- rownames(df_rownames_and_colnames)[row_index]
    # the column names are put in the second column
    col2_values <- colnames(df_rownames_and_colnames)
    # finally, we also need the actual values in this row, for the colnames
    col3_values <- as.vector(unlist(df_rownames_and_colnames[row_index, ]))
    # since our new column one is the same for the original row, we need to repeat that as many times as the columns in the original df
    col1_values <- rep(col1_value, times = length(col2_values))
    # and now put this in a new df
    slim_df_from_row <- data.frame(x = col1_values, y = col2_values, z = col3_values, stringsAsFactors = F)
    # add to larger df
    if(is.null(slim_df)){
      slim_df <- slim_df_from_row
    }
    else{
      slim_df <- rbind(slim_df, slim_df_from_row)
    }
  }
  # set the column names we have as function parameter
  colnames(slim_df) <- new_col_names
  return(slim_df)
}


get_total_interactions <- function(interactions_per_ct_list){
  # get the number of sender and receiver connections
  interaction_numbers <- get_interaction_numbers(interactions_per_ct_list)
  # get the total number of connections per cell type
  cell_types <- unique(rownames(interaction_numbers), colnames(interaction_numbers))
  # store the result in a single column dataframe
  connection_numbers <- matrix(, nrow=length(cell_types), ncol = 1, dimnames = list(cell_types, c('connections')))
  # check each cell type
  for(cell_type in cell_types){
    # nichenet does not do communication between cells of the same cell type, so we need to skip over those (NAs)
    cell_types_other <- setdiff(cell_types, cell_type)
    # get the sum of the receiving connections (rows)
    receiving_number <- sum(interaction_numbers[cell_types_other, cell_type])
    # get the sum of the sent connections (columns)
    sending_mumber <- sum(interaction_numbers[cell_type, cell_types_other])
    # sum the receiving and sending
    total_number <- receiving_number + sending_mumber
    # add this as an entry
    connection_numbers[cell_type, 'connections'] <- total_number
  }
  return(data.frame(connection_numbers))
}


get_interaction_numbers <- function(interactions_per_ct_list){
  # results are saved in a table
  interaction_numbers <- NULL
  # check the number of receptions incoming
  for(cell_type_receiver in names(interactions_per_ct_list)){
    # grab the output of this receiver
    receiver <- interactions_per_ct_list[[cell_type_receiver]]
    # and check who were the senders
    sender_cell_types <- names(receiver)
    # create a column to store the data in
    column_cell_type_receiver <- matrix(, ncol = 1, nrow = length(sender_cell_types), dimnames = list(sender_cell_types, c(cell_type_receiver)))
    # check each sender for this receiver
    for(cell_type_sender in names(receiver)){
      # extract the data for this sender
      sender <- receiver[[cell_type_sender]]
      # get the number of non-zero ligand-receptor pairs
      number_receiving_connection <- sum(sender$ligand_receptor_network > 0)
      # store this in the row that we created
      column_cell_type_receiver[cell_type_sender, cell_type_receiver] <- number_receiving_connection
    }
    # convert to dataframe for easier pasting
    column_cell_type_receiver <- data.frame(column_cell_type_receiver)
    # merge with existing data, depending on if the df was made before
    if(is.null(interaction_numbers)){
      interaction_numbers <- column_cell_type_receiver
    }
    else{
      interaction_numbers <- merge(interaction_numbers, column_cell_type_receiver, by = 0, all= T)
      # for some reason merge makes the rownames a new column, let's undo this
      rownames(interaction_numbers) <- interaction_numbers$Row.names
      interaction_numbers$Row.names <- NULL
    }
  }
  return(data.frame(interaction_numbers))
}





# locations of the files
interactions_loc <- '/data/cardiology/cell_cell_interactions/nichenet/objects/'
interactions_Baseline_t24h_loc <- paste(interactions_loc, 'v2_Baseline_vs_t24h_nichenet_onlymajors_perct.rds', sep = '')

# read object
interactions_Baseline_t24h <- readRDS(interactions_Baseline_t24h_loc)