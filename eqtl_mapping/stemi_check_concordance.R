###########################################################################################################################
#
# Libraries
#
###########################################################################################################################

require(scales)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

# plot concordance
concordance_plot <- function(x, y, ylim=c(-13,13), xlim=c(-200,200), ...) {
  # get the rows of X that are in Y as well
  x.union <- x[paste(x$SNPName, x$ProbeName) %in% paste(y$SNPName, y$ProbeName), ]
  # get the rows of Y that are in the union (match to get correct order as well)
  y.matched <- y[match(paste(x.union$SNPName, x.union$ProbeName), paste(y$SNPName, y$ProbeName)),]
  # Flip z-score if AlleleAssessed is different.
  y.matched[y.matched$AlleleAssessed != x.union$AlleleAssessed,]$OverallZScore <- y.matched[y.matched$AlleleAssessed != x.union$AlleleAssessed,]$OverallZScore * -1
  # check which are significant
  y.sign <- y.matched[y.matched$FDR < 0.05, ]
  x.sign <- x.union[x.union$FDR < 0.05, ]
  # get the union of the significant results
  x.union.overlap <- x.union[match(paste(y.sign$SNPName, y.sign$ProbeName), paste(x.sign$SNPName, x.sign$ProbeName)),]
  # calculate the concordance (same direction of significant effects)
  concordance.perc <- signif(sum(y.sign$OverallZScore * x.union.overlap$OverallZScore > 0) / length(x.union.overlap$OverallZScore) * 100, 3)
  # create empty plot
  plot.new()
  plot.window(ylim=ylim, xlim=xlim, ...)
  # create the blocks where the direction is the same, and color these green (top right and top left)
  rect(0, 2.774223, xlim[2]*2, 20, col = alpha("green", alpha = 0.2))
  rect(xlim[1]*2, -2.774223, 0, -20, col = alpha("green", alpha = 0.2))
  #create the blocks where the direction is different, and color these red (top left and bottom right)
  rect(0, -2.774223, 250, -20, col = alpha("red", alpha = 0.2))
  rect(-250, 2.774223, 0, 20, col = alpha("red", alpha = 0.2))
  # plot the points with the Z scores from the different sets being x and y
  points(x.union$OverallZScore, y.matched$OverallZScore,
         pch=19, cex=0.3, col = rgb(0,0,0, alpha = 0.2))
  # show the concordance as a legend in the plot
  text(paste0("Concordance: ", concordance.perc ,"%"), y=12, x = 120)
  # draw the rectangle in the middle to color insignificant results gray
  rect(xlim[1]*2, -2.774223, xlim[2]*2, 2.774223, col = alpha("white", alpha = 0.6), border = NA)
  # draw the lines
  abline(h=0)
  abline(v=0)
  abline(h=2.774223, lty=5)
  abline(h=-2.774223, lty=5)
  
}

###########################################################################################################################
#
# Main code
#
###########################################################################################################################

# read command line arguments
args = commandArgs(trailingOnly=TRUE)
# these we need
set1_loc <- args[1]
set2_loc <- args[2]
save_loc <- args[3]
# these are all optional
res_vert <- 600
if(!is.na(match('-res_vert', args))){
  res_vert <- args[match('-res_vert', args) + 1]
}
res_hori <- 600
if(!is.na(match('-res_hori', args))){
  res_hori <- args[match('-res_hori', args) + 1]
}
xlab = set1_loc
if(!is.na(match('-set1_name', args))){
  xlab <- args[match('-set1_name', args) + 1]
}
ylab = set2_loc
if(!is.na(match('-set2_name', args))){
  ylab <- args[match('-set2_name', args) + 1]
}
# read the two files
set1 <- read.table(set1_loc, header = T, sep = "\t", stringsAsFactors = F)
set2 <- read.table(set2_loc, header = T, sep = "\t", stringsAsFactors = F)

# set save location and resolution
png(filename = save_loc, height = res_vert, width = res_hori)
# put in the work
concordance_plot(set1, set2, xlab = xlab, ylab = ylab)
dev.off()


