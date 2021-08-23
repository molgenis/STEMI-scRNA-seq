

# functions

eqtl_result_to_plottable_table <- function(eqtl_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'UT_Baseline', 'UT_t24h', 'UT_t8w')){
  # init the table
  plot_table <- NULL
  # check each condition
  for(condition in conditions){
    # check each cell type
    for(cell_type in cell_types){
      # paste full location together
      full_output_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/', 'eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      # set default
      nr_of_egenes <- 0
      try({
        # read the file
        eqtl_output <- read.table(full_output_loc, sep = '\t', header = T, stringsAsFactors = F)
        # get the number of egenes
        nr_of_egenes <- length(unique(eqtl_output[eqtl_output$FDR < 0.05, 'ProbeName']))
      })
      # turn into a row
      this_row <- data.frame(condition=c(condition), cell_type=c(cell_type), egenes=c(nr_of_egenes), stringsAsFactors = F)
      # add to table
      if(is.null(plot_table)){
        plot_table <- this_row
      }
      else{
        plot_table <- rbind(plot_table, this_row)
      }
    }
  }
  return(plot_table)
}

plot_eqtl_numbers_per_condition <- function(eqtl_table, use_label_dict=T, use_color_dict=T, paper_style=F){
  # get pretty labels if possible
  if(use_label_dict){
    labels_dict <- label_dict()
    eqtl_table$cell_type <- as.vector(unlist(labels_dict[as.character(eqtl_table$cell_type)]))
    eqtl_table$condition <- as.vector(unlist(labels_dict[as.character(eqtl_table$condition)]))
  }
  # we want everything to be of the same size
  ylim <- max(eqtl_table$egenes)
  # make a list of plots
  plot_list <- list()
  # store the legend
  legend <- NULL
  # make a plot per condition
  for(condition in unique(eqtl_table$condition)){
    print(head(eqtl_table[eqtl_table$condition == condition, ]))
    # make the plot
    p <- ggplot(data=eqtl_table[eqtl_table$condition == condition, ], aes(x=cell_type, y=egenes, fill=cell_type)) +
      geom_bar(stat='identity') +
      ggtitle(paste('number of egenes in', condition)) +
      xlab('cell type') +
      ylab('number of egenes')
    # make a color scale
    if(use_color_dict){
      # fetch colors
      cc <- get_color_coding_dict()
      # set colors based on condition
      colScale <- scale_fill_manual(name = "cell type",values = unlist(cc[eqtl_table[eqtl_table$condition == condition, 'cell_type']]))
      # add color
      p <- p + colScale
    }
    # set the y limit
    p <- p + ylim(c(0, ylim))
    # extract legend
    legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12)))
    # remove the legend from the plot
    p <- p + theme(legend.position = 'none')
    # apply the minimalist paper style (bleh)
    if(paper_style){
      p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    }
    # add to plot list
    plot_list[[condition]] <- p
  }
  # add one legend
  #plot_list[['legend']] <- legend
  # plot in grid
  final_plot <- plot_grid(plotlist = plot_list)
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
  color_coding[["HC-t24h"]] <- "khaki4"
  color_coding[["HC-t8w"]] <- "paleturquoise1"
  color_coding[["t0-t24h"]] <- "paleturquoise3"
  color_coding[["t0-t8w"]] <- "rosybrown1"
  color_coding[["t24h-t8w"]] <- "rosybrown3"
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
  # eqtl combined combinations
  label_dict[['UT_Baseline']] <- 'HC+t0'
  label_dict[['UT_t24h']] <- 'HC+t24h'
  label_dict[['UT_t8w']] <- 'HC+t8w'
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
  return(label_dict)
}


# locations
eqtl_runs_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/inhouse_eQTL_mapping_pipeline/'
eqtl_output_loc <- paste(eqtl_runs_loc, 'stemi_and_1mut_meta_lowerres_20210629/', sep = '')
eqtl_output_loc_eqtlgenlead <- paste(eqtl_runs_loc, 'stemi_and_1mut_meta_lowerres_20210629_confine_lead_snp_gene/', sep = '')

# get plottable table
eqtl_number_table <- eqtl_result_to_plottable_table(eqtl_output_loc)
# make the plot
plot_eqtl_numbers_per_condition(eqtl_number_table)
# get plottable table
eqtl_number_table_eqtlgend_lead <- eqtl_result_to_plottable_table(eqtl_output_loc_eqtlgenlead)
# make the plot
plot_eqtl_numbers_per_condition(eqtl_number_table_eqtlgend_lead)


