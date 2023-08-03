#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_plot_lr_connection_circles.R
# Function: 
############################################################################################################################

####################
# libraries        #
####################
library(circlize)


####################
# Functions        #
####################

interactions_to_circle <- function(interactions_per_ct_list, plot_title='cell communication', by_receptor=T, cutoff=0.05, split_sender_and_receiver=F, order_alphabet=T, order_size=F, start_degree=90){
  # get the total number of connections per cell type
  interaction_numbers <- get_interaction_numbers(interactions_per_ct_list, by_receptor, cutoff)
  # slim it down into a three column input
  connections_slim <- slim_df_down(interaction_numbers, new_col_names = c('sender', 'receiver', 'connections'))
  # split the sender and receiver is requested
  if (split_sender_and_receiver) {
    connections_slim[['sender']] <- paste(connections_slim[['sender']], 'sent')
    connections_slim[['receiver']] <- paste(connections_slim[['receiver']], 'received')
  }
  # remove empty entries
  connections_slim <- connections_slim[!is.na(connections_slim[['connections']]) & connections_slim[['connections']] > 0, ]
  # add receiver and sender together
  connections_all_ct <- combine_directions_cells(connections_slim)
  # sort on size
  if (order_alphabet) {
    connections_all_ct <- connections_all_ct[order(connections_all_ct[['cell_type']]), ]
  }
  if (order_size) {
    connections_all_ct <- connections_all_ct[order(connections_all_ct[['connections']]), ]
  }
  # is we split the sender and receiver, we need to split those up
  if (split_sender_and_receiver) {
    connections_all_ct <- rbind(connections_all_ct[grepl(' received$', connections_all_ct[['cell_type']]), ],
          connections_all_ct[grepl(' sent$', connections_all_ct[['cell_type']]), ])
  }
  # okay, I did this somewhere in the past, so afraid I have to do it again
  colnames(connections_all_ct) <- c('a', 'b')
  # get the data in a start/stop manner, instead of a size-manner
  connections_start_stop <- turn_sizes_to_ranges(connections_slim, 'connections', 'sender', 'receiver')
  # make sure the other plot is gone
  circos.clear()
  # set the start degree
  circos.par(start.degree=start_degree)
  # start making the plot, set height of track (cell types)
  circos.par('track.height' = 0.1)
  xlims <- data.frame(a=rep(0, times=nrow(connections_all_ct)), b=connections_all_ct$b)
  # start building
  circos.initialize(sectors = connections_all_ct$a, xlim = xlims)
  # add track with labels
  circos.track(connections_all_ct$a, y=connections_all_ct$b, panel.fun = function(x, y){
    circos.text(CELL_META$xcenter, 
                CELL_META$ycenter, 
                CELL_META$sector.index)
    circos.axis(labels.cex = 0.6, labels.facing = 'outside')
  })
  # color the tracks
  for(cell_type in unique(connections_all_ct$a)){
    draw.sector(get.cell.meta.data("cell.start.degree", sector.index = cell_type),
                get.cell.meta.data("cell.end.degree", sector.index = cell_type),
                rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
                rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1),
                col = get_color_coding_dict()[[cell_type]])
  }
  # draw colour on the sectors, representing the cell types
  for(cell_type in unique(connections_all_ct$a)){
    highlight.sector(c(cell_type), track.index = 1, text=label_dict()[[cell_type]], col = '#ffffff00', text.col = '#ffffffff')
  }
  # draw the connections
  apply(connections_start_stop, 1, function(row){
    if(!is.na(row['sender'])){
      # grab the color of the sender
      send_color <- get_color_coding_dict()[[row['sender']]]
      # and of the receiver
      receiver_color <- get_color_coding_dict()[[row['receiver']]]
      # split colors in five
      ramp_colors <- colorRampPalette(c(send_color, receiver_color))(5)
      # color slightly more towards the sender, by taking the second color from the 5 levels, and add remove possible existing transparancy
      connection_color <- substr(ramp_colors[2], 1, 7)
      # manually add transparancy by adding it to normal 6-colour hex code
      connection_color <- paste(connection_color, '80', sep = '')
      # grab the positions on the sender
      from_start <- as.numeric(row['from_start'])
      from_stop <- as.numeric(row['from_stop'])
      # grab the positions on the receiver
      to_start <- as.numeric(row['to_start'])
      to_stop <- as.numeric(row['to_stop'])
      # draw the connection
      circos.link(row['sender'], c(from_start, from_stop), row['receiver'], c(to_start, to_stop), col = connection_color, directional = 1)
    }
  })
  title(plot_title)
  circos.clear()
}


interactions_to_percentages <- function(interaction_numbers, percentage_of='both_receiver'){
  # create a new matrix with the same dimensions as the original
  interaction_percentages <- matrix(, ncol=ncol(interaction_numbers), nrow=nrow(interaction_numbers), dimnames = list(rownames(interaction_numbers), colnames(interaction_numbers)))
  # calculate the sum of the rows, which are the senders
  sent_sums <- as.list(rowSums(interaction_numbers, na.rm = T))
  names(sent_sums) <- rownames(interaction_numbers)
  # calculate the sum of the columns, which are the receivers
  received_sums <- as.list(colSums(interaction_numbers, na.rm = T))
  names(received_sums) <- colnames(interaction_numbers)
  # check each sender
  for(sender in rownames(interaction_numbers)){
    # check each receiver
    for(receiver in colnames(interaction_numbers)){
      # grab the value for that combination
      sender_receiver_combo  <- interaction_numbers[sender, receiver]
      # depending on the choice, we will calculate a percentage of senders, receivers, or all
      if(percentage_of == 'sender'){
        sender_receiver_combo <- sender_receiver_combo / sent_sums[[receiver]]
      }
      else if(percentage_of == 'receiver'){
        sender_receiver_combo <- sender_receiver_combo / received_sums[[receiver]]
      }
      else if(percentage_of == 'both_sender'){
        sender_receiver_combo <- sender_receiver_combo / (sent_sums[[sender]] + received_sums[[sender]])
      }
      else if(percentage_of == 'both_receiver'){
        sender_receiver_combo <- sender_receiver_combo / (sent_sums[[receiver]] + received_sums[[receiver]])
      }
      else{
        print('unknown option, doing both receiver')
        sender_receiver_combo <- sender_receiver_combo / (sent_sums[[receiver]] + received_sums[[receiver]])
      }
      # set the value
      interaction_percentages[sender, receiver] <- sender_receiver_combo
    }
  }
  return(data.frame(interaction_percentages))
}


all_interactions_to_percentages <- function(interactions_per_ct_list, output_loc_prepend, by_receptor=T, percentages_of=c('sender', 'receiver', 'both_sender', 'both_receiver')){
  # first get the numbers
  interaction_numbers <- get_interaction_numbers(interactions_per_ct_list, by_receptor)
  # then calculate each percentage
  for(percentage_of in percentages_of){
    interaction_percentage <- interactions_to_percentages(interaction_numbers, percentage_of)
    # paste together an output path
    full_out <- paste(output_loc_prepend, '_percentage_', percentage_of, '.tsv', sep = '')
    # write the result
    write.table(interaction_percentage, full_out, sep = '\t', quote = F, col.names = T)
  }
}


remove_empty_entries <- function(sender_receiver_table){
  # get the possible cell types
  cell_types <- unique(rownames(sender_receiver_table), colnames(sender_receiver_table))
  # initialize values
  n_ct_connections <- as.list(rep(0, length(cell_types)))
  names(n_ct_connections) <- cell_types
  # check each entry
  for(ct_a in rownames(sender_receiver_table)){
    for(ct_b in colnames(sender_receiver_table)){
      # check for this combination
      combination <- sender_receiver_table[ct_a, ct_b]
      # if it is either a sender or a receiver, they have a connection
      if(!is.na(combination) & combination > 0){
        n_ct_connections[[ct_a]] <-  n_ct_connections[[ct_a]] + combination
        n_ct_connections[[ct_b]] <-  n_ct_connections[[ct_b]] + combination
      }
    }
  }
  # grab the entries that were not zero
  non_zeroes <- names(n_ct_connections[n_ct_connections != 0])
  # subset the table
  subsetted_table <- sender_receiver_table[non_zeroes, non_zeroes]
  return(subsetted_table)
}



turn_sizes_to_ranges <- function(link_sizes_df, size_column, from_column, to_column){
  # get the possible entries
  entries <- unique(c(link_sizes_df[[from_column]], link_sizes_df[[to_column]]))
  # create the dataframe to store everything in
  link_df <- data.frame(matrix(, nrow=(length(entries))*(length(entries)-1)/2, ncol=6))
  colnames(link_df) <- c(from_column, to_column, 'from_start', 'from_stop','to_start', 'to_stop')
  # we'll add entries by index, that is a lot faster than rbinding new rows
  i <- 1
  # check each entry
  for(x in 1:(length(entries))){
    # the column on the horizontal axis is the number, and we'll grab the entry for that number
    from_entry <- entries[x]
    # from the bottom left and from the top right are the same, and the entry against itself is always empty
    for(y in (1):(length(entries))){
      # grab on the vertical axis
      to_entry <- entries[y]
      # store the link size
      link_size <- NA
      # get the entry if it exists
      if(nrow(link_sizes_df[link_sizes_df[[from_column]] == from_entry & link_sizes_df[[to_column]] == to_entry, ]) > 0){
        link_size <- as.numeric(link_sizes_df[link_sizes_df[[from_column]] == from_entry & link_sizes_df[[to_column]] == to_entry, size_column])
      }
      # if this link exists, we can continue
      if(!is.na(link_size)){
        max_from <- max(c(0,
                          as.numeric(link_df[!is.na(link_df[[from_column]]) & !is.na(link_df[[to_column]]) &
                                               link_df[[from_column]] == from_entry, 'from_stop']),
                          as.numeric(link_df[!is.na(link_df[[from_column]]) & !is.na(link_df[[to_column]]) &
                                               link_df[[to_column]] == from_entry, 'to_stop'])
        )
        )
        max_to <- max(c(0,
                        as.numeric(link_df[!is.na(link_df[[to_column]]) & !is.na(link_df[[to_column]]) &
                                             link_df[[from_column]] == to_entry, 'from_stop']),
                        as.numeric(link_df[!is.na(link_df[[to_column]]) & !is.na(link_df[[to_column]]) &
                                             (link_df[[to_column]] == to_entry), 'to_stop'])
        )
        )
        
        # now we are either zero or the max of both entries, we'll turn that into the actual entries
        link_df[i, from_column] <- from_entry
        link_df[i, to_column] <- to_entry
        link_df[i, 'from_start'] <- max_from
        link_df[i, 'from_stop'] <- (max_from + link_size)
        link_df[i, 'to_start'] <- max_to
        link_df[i, 'to_stop'] <- (max_to + link_size)
        # update the index
        i <- i + 1
      }
    }
  }
  return(link_df)
}


combine_directions_cells <- function(slim_connection_df, cell_types){
  # we will save the sum of the incoming and outgoing connections per cell type
  total_sum_per_cell_type <- NULL
  # we will check each sender
  for(cell_type in unique(c(slim_connection_df[['sender']], slim_connection_df[['receiver']]))){
    # check that sender against each receiver
    sent_connections <- sum(slim_connection_df[!is.na(slim_connection_df[['connections']]) & slim_connection_df[['sender']] == cell_type, 'connections'])
    # check the incoming the connections
    incoming_connections <- sum(slim_connection_df[!is.na(slim_connection_df[['connections']]) & slim_connection_df[['receiver']] == cell_type, 'connections'])
    # add these two together
    total_connections <- sent_connections + incoming_connections
    # turn this into a row for tht dataframe
    connection_row <- data.frame(cell_type=c(cell_type), connections=total_connections)
    # add to dataframe
    if(is.null(total_sum_per_cell_type)){
      total_sum_per_cell_type <- connection_row
    }
    else{
      total_sum_per_cell_type <- rbind(total_sum_per_cell_type, connection_row)
    }
  }
  return(total_sum_per_cell_type)
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


get_interaction_numbers <- function(interactions_per_ct_list, by_receptor=T, cutoff=0.05){
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
      # init variable
      number_receiving_connection <- NA
      # check if the user want to use the targets or the receptors as connections to the ligand
      if(by_receptor){
        # get the number of non-zero ligand-receptor pairs
        number_receiving_connection <- sum(sender$ligand_receptor_network > cutoff)
      }
      else{
        # get the number of non-zero ligand-target pairs
        number_receiving_connection <- sum(sender$ligand_target > cutoff)
      }
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


plots_to_files <- function(plots_per_pathway, path_prepend='./', title_prepend='interactions', ncols=2){
  # paste together the path to the output
  output_path <- paste(path_prepend, 'interactions.pdf', sep = '')
  # check how many pathways there are
  n_pathways <- length(names(plots_per_pathway))
  n_row_frac <- n_pathways/ncols
  n_row <- ceiling(n_row_frac)
  # set up the  rows and columns
  par(mfrow=c(n_row,ncols))
  # we'll save this output
  pdf(output_path)
  
  # do per pathway
  for(pathway in names(plots_per_pathway)){
    # get this specific pathway
    plots_pathway <- plots_per_pathway[[pathway]]
    tryCatch({
      # turn into plot
      interactions_to_circle(plots_pathway, plot_title = paste(title_prepend, pathway), by_receptor = F)
    }, error=function(cond) {
      print(paste('got error for pathways', pathway))
      print(cond)
    }
    )
    
  }
  # these were all the plots
  dev.off()
}


get_ligand_target_combinations <- function(nichenet_object){
  # get the ligands and the downstream targets
  ligands <- rownames(nichenet_object$ligand_target)
  targets <- colnames(nichenet_object$ligand_target)
  # store the interacting components
  interactions <- c()
  # check each ligand
  for(ligand in ligands){
    # check each target
    for(target in targets){
      interaction <- nichenet_object$ligand_target[ligand, target]
      # if there is an interaction
      if(interaction > 0){
        # add it to our interactions
        interactions <- c(interactions, paste(ligand, target, sep = '_'))
      }
    }
  }
  return(interactions)
}


check_concordance_nichenet_objects <- function(nichenet_object_1, nichenet_object_2){
  # get the interactions from the first object
  interactions_1 <- get_ligand_target_combinations(nichenet_object_1)
  # get the interactions from the second object
  interactions_2 <- get_ligand_target_combinations(nichenet_object_2)
  # get the number that is shared
  shared_size <- length(intersect(interactions_1, interactions_2))
  # get what is unique in 1
  unique_1 <- length(setdiff(interactions_1, interactions_2))
  # get what is unique in 2
  unique_2 <- length(setdiff(interactions_2, interactions_1))
  # return result
  return(list('shared'=shared_size, 'unique_1'=unique_1, 'unique_2'=unique_2))
}


check_concordance_nichenet_objects_per_ct_combination <- function(nichenet_list1, nichenet_list2){
  # put the result in a table
  sharing_table <- NULL
  # check which cell types we can do
  cell_types_senders_to_check <- intersect(names(nichenet_list1), names(nichenet_list2))
  # check each of these cell types
  for(cell_type_sending in cell_types_senders_to_check){
    # subset to that
    nichenet_receivers_list_1 <- nichenet_list1[[cell_type_sending]]
    nichenet_receivers_list_2 <- nichenet_list2[[cell_type_sending]]
    # now check which of these overlap
    cell_types_receivers_to_check <- intersect(names(nichenet_receivers_list_1), names(nichenet_receivers_list_2))
    # check each of these cell types
    for(cell_type_receiving in cell_types_receivers_to_check){
      # check for a result
      if(length(nichenet_receivers_list_1[[cell_type_receiving]])>0 & length(nichenet_receivers_list_2[[cell_type_receiving]])>0){
        # get the sharing
        sharing <- check_concordance_nichenet_objects(nichenet_receivers_list_1[[cell_type_receiving]], nichenet_receivers_list_2[[cell_type_receiving]])
        # turn list into dataframe
        sharing <- as.data.frame(sharing)
        # add receiving cell type
        sharing[['receiver']] <- cell_type_receiving
        sharing[['sender']] <- cell_type_sending
        # add to the table
        if(is.null(sharing_table)){
          sharing_table <- sharing
        }
        else{
          sharing_table <- rbind(sharing_table, sharing)
        }
      }
    }
    # add the sending cell type
  }
  return(sharing_table)
}

check_concordance_nichenet_objects_per_pathway <- function(nichenet_list1, nichenet_list2){
  # we will save the result
  result_table <- NULL
  # check which pathways we can use
  common_pathways <- intersect(names(nichenet_list1), names(nichenet_list2))
  # check these
  for(pathway in common_pathways){
    # grab for that specific pathway
    nichenet_list_pathway_1 <- nichenet_list1[[pathway]]
    nichenet_list_pathway_2 <- nichenet_list2[[pathway]]
    # check for this specific pathway
    sharing <- check_concordance_nichenet_objects_per_ct_combination(nichenet_list_pathway_1, nichenet_list_pathway_2)
    # add the pathway name
    sharing[['pathway']] <- pathway
    # add to result table
    if(is.null(result_table)){
      result_table <- sharing
    }
    else{
      result_table <- rbind(result_table, sharing)
    }
  }
  return(result_table)
}


get_interactions_from_table <- function(interaction_table, cutoff=0.05){
  # we have a vector in which we will store the interactions
  interactions <- c()
  # check each of the downstream genes
  downstream_genes <- rownames(interaction_table)
  # check each of the ligands
  ligands <- colnames(interaction_table)
  # now check each combination
  for(downstream_gene in downstream_genes){
    for(ligand in ligands){
      # check the interaction
      interaction <- interaction_table[downstream_gene, ligand]
      # check if the interaction is large enough
      if(interaction > cutoff){
        # add it to the list of interactions
        interactions <- c(interactions, paste(ligand, downstream_gene, sep = '_'))
      }
    }
  }
  return(interactions)
}

combine_interactions_per_pathways <- function(interactions_per_pathway_1, interactions_per_pathway_2, cutoff=0.05, by_receptor=T){
  # store the interactions per pathway
  interactions_per_pathway_combined <- list()
  # check which pathways we can use
  common_pathways <- intersect(names(interactions_per_pathway_1), names(interactions_per_pathway_2))
  # check these
  for(pathway in common_pathways){
    # we will save the result
    interaction_numbers <- NULL
    # grab for that specific pathway
    nichenet_list_pathway_1 <- interactions_per_pathway_1[[pathway]]
    nichenet_list_pathway_2 <- interactions_per_pathway_2[[pathway]]
    # check which cell types we can do
    cell_types_senders_to_check <- intersect(names(nichenet_list_pathway_1), names(nichenet_list_pathway_2))
    # check each of these cell types
    for(cell_type_sending in cell_types_senders_to_check){
      # subset to that
      nichenet_receivers_list_1 <- nichenet_list_pathway_1[[cell_type_sending]]
      nichenet_receivers_list_2 <- nichenet_list_pathway_2[[cell_type_sending]]
      # now check which of these overlap
      cell_types_receivers_to_check <- intersect(names(nichenet_receivers_list_1), names(nichenet_receivers_list_2))
      # check each of these cell types
      for(cell_type_receiving in cell_types_receivers_to_check){
        # grab the interactions
        interactions_1 <- NULL
        interactions_2 <- NULL
        if(by_receptor){
          interactions_1 <- nichenet_receivers_list_1[[cell_type_receiving]][['ligand_receptor_network']]
          interactions_2 <- nichenet_receivers_list_2[[cell_type_receiving]][['ligand_receptor_network']]
        }
        else{
          interactions_1 <- nichenet_receivers_list_1[[cell_type_receiving]][['ligand_target']]
          interactions_2 <- nichenet_receivers_list_2[[cell_type_receiving]][['ligand_target']]
        }
        # get the interactions from the tables
        interactions_genes_1 <- get_interactions_from_table(interactions_1, cutoff = cutoff)
        interactions_genes_2 <- get_interactions_from_table(interactions_2, cutoff = cutoff)
        # get the outer join of the genes
        interactions_genes <- unique(c(interactions_genes_1, interactions_genes_2))
        # count the number of interactions
        interaction_number <- length(interactions_genes)
        # create a row
        interaction_row <- data.frame(from=c(cell_type_sending), to=c(cell_type_receiving), number=c(interaction_number))
        # if this is the first round through, we will initialize the dataframe
        if(is.null(interaction_numbers)){
          interaction_numbers <- interaction_row
        }
        else{
          interaction_numbers <- rbind(interaction_numbers, interaction_row)
        }
      }
    }
    interactions_per_pathway_combined[[pathway]] <- interaction_numbers
  }
  return(interactions_per_pathway_combined)
}


interactions_to_circle_from_table <- function(interactions_table, plot_title='cell communication'){
  # set the correct column names
  colnames(interactions_table) <- c('sender', 'receiver', 'connections')
  # remove the empty connections
  interactions_table <- interactions_table[interactions_table[['connections']] > 0, ]
  # add receiver and sender together
  connections_all_ct <- combine_directions_cells(interactions_table)
  # okay, I did this somewhere in the past, so afraid I have to do it again
  colnames(connections_all_ct) <- c('a', 'b')
  # get the data in a start/stop manner, instead of a size-manner
  connections_start_stop <- turn_sizes_to_ranges(interactions_table, 'connections', 'sender', 'receiver')
  # make sure the other plot is gone
  circos.clear()
  # start making the plot, set height of track (cell types)
  circos.par('track.height' = 0.1)
  xlims <- data.frame(a=rep(0, times=nrow(connections_all_ct)), b=connections_all_ct$b)
  # start building
  circos.initialize(sectors = connections_all_ct$a, xlim = xlims)
  # add track with labels
  circos.track(connections_all_ct$a, y=connections_all_ct$b, panel.fun = function(x, y){
    circos.text(CELL_META$xcenter, 
                CELL_META$ycenter, 
                CELL_META$sector.index)
    circos.axis(labels.cex = 0.6, labels.facing = 'outside')
  })
  # color the tracks
  for(cell_type in unique(connections_all_ct$a)){
    draw.sector(get.cell.meta.data("cell.start.degree", sector.index = cell_type),
                get.cell.meta.data("cell.end.degree", sector.index = cell_type),
                rou1 = get.cell.meta.data("cell.top.radius", track.index = 1),
                rou2 = get.cell.meta.data("cell.bottom.radius", track.index = 1),
                col = get_color_coding_dict()[[cell_type]])
  }
  # draw colour on the sectors, representing the cell types
  for(cell_type in unique(connections_all_ct$a)){
    highlight.sector(c(cell_type), track.index = 1, text=label_dict()[[cell_type]], col = '#ffffff00', text.col = '#ffffffff')
  }
  # draw the connections
  apply(connections_start_stop, 1, function(row){
    if(!is.na(row['sender'])){
      # grab the color of the sender
      send_color <- get_color_coding_dict()[[row['sender']]]
      # and of the receiver
      receiver_color <- get_color_coding_dict()[[row['receiver']]]
      # split colors in five
      ramp_colors <- colorRampPalette(c(send_color, receiver_color))(5)
      # color slightly more towards the sender, by taking the second color from the 5 levels, and add remove possible existing transparancy
      connection_color <- substr(ramp_colors[2], 1, 7)
      # manually add transparancy by adding it to normal 6-colour hex code
      connection_color <- paste(connection_color, '80', sep = '')
      # grab the positions on the sender
      from_start <- as.numeric(row['from_start'])
      from_stop <- as.numeric(row['from_stop'])
      # grab the positions on the receiver
      to_start <- as.numeric(row['to_start'])
      to_stop <- as.numeric(row['to_stop'])
      # draw the connection
      circos.link(row['sender'], c(from_start, from_stop), row['receiver'], c(to_start, to_stop), col = connection_color, directional = 1)
    }
  })
  title(plot_title)
}


plots_to_files_from_tables <- function(communication_tables_per_pathways, path_prepend='./', title_prepend='interactions', ncols=2){
  # paste together the path to the output
  output_path <- paste(path_prepend, 'interactions.pdf', sep = '')
  # check how many pathways there are
  n_pathways <- length(names(communication_tables_per_pathways))
  n_row_frac <- n_pathways/ncols
  n_row <- ceiling(n_row_frac)
  # set up the  rows and columns
  par(mfrow=c(n_row,ncols))
  # we'll save this output
  pdf(output_path)
  # check each pathway
  for(pathway in names(communication_tables_per_pathways)){
    tryCatch({
      # get the specific table
      communications_table <- communication_tables_per_pathways[[pathway]]
      # create the plot
      interactions_to_circle_from_table(communications_table, pathway)
    }, error=function(cond) {
      print(paste('got error for pathways', pathway))
      print(cond)
    }
    )
    
  }
  # these were all the plots
  dev.off()
}


tables_to_files <- function(communication_numbers_per_pathway, prepend){
  # check each pathway
  for(pathway in names(communication_numbers_per_pathway)){
    # set up the output location
    output_loc <- paste(prepend, pathway, '.tsv', sep = '')
    # write the result
    write.table(communication_numbers_per_pathway[[pathway]], output_loc, sep = '\t', col.names = T, row.names = F)
  }
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
  color_coding[["UT\nBaseline"]] <- "khaki2"
  color_coding[["UT\nt24h"]] <- "khaki4"
  color_coding[["UT\nt8w"]] <- "paleturquoise1"
  color_coding[["Baseline\nt24h"]] <- "paleturquoise3"
  color_coding[["Baseline\nt8w"]] <- "rosybrown1"
  color_coding[["t24h\nt8w"]] <- "rosybrown3"
  color_coding[["UT-Baseline"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["Baseline-t24h"]] <- "paleturquoise3"
  color_coding[["Baseline-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
  color_coding[["UT-t0"]] <- "khaki2"
  color_coding[["UT-t24h"]] <- "khaki4"
  color_coding[["UT-t8w"]] <- "paleturquoise1"
  color_coding[["HC-t0"]] <- "khaki2"
  color_coding[["t0-HC"]] <- "khaki2"
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["t24h-HC"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t8w-HC"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t24h-t0"]] <- "#FF6066" #"paleturquoise3"
  color_coding[["t0-t8w"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t8w-t0"]] <- "#C060A6" #"rosybrown1"
  color_coding[["t24h-t8w"]] <- "#C00040" #"rosybrown3"
  color_coding[["t8w-t24h"]] <- "#C00040" #"rosybrown3"
  # set condition colors
  color_coding[["HC"]] <- "grey"
  color_coding[["t0"]] <- "pink"
  color_coding[["t24h"]] <- "red"
  color_coding[["t8w"]] <- "purple"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # sent and received split
  color_coding[["CD4T sent"]] <- "#153057"
  color_coding[["CD8T sent"]] <- "#009DDB"
  color_coding[["monocyte sent"]] <- "#EDBA1B"
  color_coding[["NK sent"]] <- "#E64B50"
  color_coding[["B sent"]] <- "#71BC4B"
  color_coding[["DC sent"]] <- "#965EC8"
  color_coding[["CD4+ T sent"]] <- "#153057"
  color_coding[["CD8+ T sent"]] <- "#009DDB"
  color_coding[["CD4T received"]] <- "#153057"
  color_coding[["CD8T received"]] <- "#009DDB"
  color_coding[["monocyte received"]] <- "#EDBA1B"
  color_coding[["NK received"]] <- "#E64B50"
  color_coding[["B received"]] <- "#71BC4B"
  color_coding[["DC received"]] <- "#965EC8"
  color_coding[["CD4+ T received"]] <- "#153057"
  color_coding[["CD8+ T received"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  # condition combinations
  label_dict[['UTBaseline']] <- 'UT-Baseline'
  label_dict[['UTt24h']] <- 'UT-t24h'
  label_dict[['UTt8w']] <- 'UT-t8w'
  label_dict[['Baselinet24h']] <- 'Baseline-t24h'
  label_dict[['Baselinet8w']] <- 'Baseline-t8w'
  label_dict[['t24ht8w']] <- 't24h-t8w'
  label_dict[['UTBaseline']] <- 'HC-t0'
  label_dict[['UTt24h']] <- 'HC-t24h'
  label_dict[['UTt8w']] <- 'HC-t8w'
  label_dict[['Baselinet24h']] <- 't0-t24h'
  label_dict[['Baselinet8w']] <- 't0-t8w'
  # conditions
  label_dict[['UT']] <- 'HC'
  label_dict[['Baseline']] <- 't0'
  label_dict[['t24h']] <- 't24h'
  label_dict[['t8w']] <- 't8w'
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  # and sender receiver split
  label_dict[["CD4T sent"]] <- "CD4+ T sent"
  label_dict[["CD8T sent"]] <- "CD8+ T sent"
  label_dict[["monocyte sent"]] <- "monocyte sent"
  label_dict[["NK sent"]] <- "NK sent"
  label_dict[["B sent"]] <- "B sent"
  label_dict[["DC sent"]] <- "DC sent"
  label_dict[["HSPC sent"]] <- "HSPC sent"
  label_dict[["plasmablast sent"]] <- "plasmablast sent"
  label_dict[["platelet sent"]] <- "platelet sent"
  label_dict[["T_other sent"]] <- "other T sent"
  label_dict[["CD4T received"]] <- "CD4+ T received"
  label_dict[["CD8T received"]] <- "CD8+ T received"
  label_dict[["monocyte received"]] <- "monocyte received"
  label_dict[["NK received"]] <- "NK received"
  label_dict[["B received"]] <- "B received"
  label_dict[["DC received"]] <- "DC received"
  label_dict[["HSPC received"]] <- "HSPC received"
  label_dict[["plasmablast received"]] <- "plasmablast received"
  label_dict[["platelet received"]] <- "platelet received"
  label_dict[["T_other received"]] <- "other T received"
  return(label_dict)
}

text_color_dict <- function(){
  text_color_dict <- list()
  # set the cell type colors
  text_color_dict[["Bulk"]] <- "#000000ff"
  text_color_dict[["CD4T"]] <- "#000000ff"
  text_color_dict[["CD8T"]] <- "#000000ff"
  text_color_dict[["monocyte"]] <- "#ffffffff"
  text_color_dict[["NK"]] <- "#000000ff"
  text_color_dict[["B"]] <- "#fffffff"
  text_color_dict[["DC"]] <- "#000000ff"
  text_color_dict[["CD4+ T"]] <- "#000000ff"
  text_color_dict[["CD8+ T"]] <- "#000000ff"
  # other cell type colors
  text_color_dict[["HSPC"]] <- "#000000ff"
  text_color_dict[["platelet"]] <- "#000000ff"
  text_color_dict[["plasmablast"]] <- "#000000ff"
  text_color_dict[["other T"]] <- "#000000ff"
  # sent and received split
  text_color_dict[["CD4T sent"]] <- "#000000ff"
  text_color_dict[["CD8T sent"]] <- "#000000ff"
  text_color_dict[["monocyte sent"]] <- "#ffffffff"
  text_color_dict[["NK sent"]] <- "#000000ff"
  text_color_dict[["B sent"]] <- "#fffffff"
  text_color_dict[["DC sent"]] <- "#000000ff"
  text_color_dict[["CD4+ T sent"]] <- "#000000ff"
  text_color_dict[["CD8+ T sent"]] <- "#000000ff"
  text_color_dict[["HSPC sent"]] <- "#000000ff"
  text_color_dict[["platelet sent"]] <- "#000000ff"
  text_color_dict[["plasmablast sent"]] <- "#000000ff"
  text_color_dict[["other T sent"]] <- "#000000ff"
  text_color_dict[["CD4T received"]] <- "#000000ff"
  text_color_dict[["CD8T received"]] <- "#000000ff"
  text_color_dict[["monocyte received"]] <- "#ffffffff"
  text_color_dict[["NK received"]] <- "#000000ff"
  text_color_dict[["B received"]] <- "#fffffff"
  text_color_dict[["DC received"]] <- "#000000ff"
  text_color_dict[["CD4+ T received"]] <- "#000000ff"
  text_color_dict[["CD8+ T received"]] <- "#000000ff"
  text_color_dict[["HSPC received"]] <- "#000000ff"
  text_color_dict[["platelet received"]] <- "#000000ff"
  text_color_dict[["plasmablast received"]] <- "#000000ff"
  text_color_dict[["other T received"]] <- "#000000ff"
  return(text_color_dict)
}


####################
# Main Code        #
####################

# locations of the files
interactions_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/cell_cell_interactions/nichenet/objects/'
interactions_Baseline_t24h_v2_loc <- paste(interactions_loc, 'v2_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted.rds', sep = '')
interactions_t24h_t8w_v2_loc <- paste(interactions_loc, 'v2_t24h_vs_t8w_nichenet_onlymajor_perct_omni_unweighted.rds', sep = '')
interactions_Baseline_t24h_v3_loc <- paste(interactions_loc, 'v3_Baseline_vs_t24h_nichenet_onlymajors_perct_omni_unweighted.rds', sep = '')
interactions_t24h_t8w_v3_loc <- paste(interactions_loc, 'v3_t24h_vs_t8w_nichenet_onlymajor_perct_omni_unweighted.rds', sep = '')


# read object
interactions_Baseline_t24h_v2 <- readRDS(interactions_Baseline_t24h_v2_loc)
interactions_t24h_t8w_v2 <- readRDS(interactions_t24h_t8w_v2_loc)
interactions_Baseline_t24h_v3 <- readRDS(interactions_Baseline_t24h_v3_loc)
interactions_t24h_t8w_v3 <- readRDS(interactions_t24h_t8w_v3_loc)

# setup four plot tiles
par(mfrow=c(2,2))
interactions_to_circle(interactions_Baseline_t24h_v2, plot_title = title('V2 Cell communication changes \nbetween t0 and t24h'), by_receptor = T, split_sender_and_receiver = T)
interactions_to_circle(interactions_Baseline_t24h_v3, plot_title = title('V3 Cell communication changes \nbetween t0 and t24h'), by_receptor = T, split_sender_and_receiver = T)
interactions_to_circle(interactions_t24h_t8w_v2, plot_title = title('V2 Cell communication changes \nbetween t24h and t8w'), by_receptor = T, split_sender_and_receiver = T)
interactions_to_circle(interactions_t24h_t8w_v3, plot_title = title('V3 Cell communication changes \nbetween t24h and t8w'), by_receptor = T, split_sender_and_receiver = T)
par(mfrow=c(1,1))
