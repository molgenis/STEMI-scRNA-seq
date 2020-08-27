library(data.table)

# convert the genotypes to numbers for correlation
convert_genotypes <- function(genotype_file){
  # convert for each genotype column
  for(gt in colnames(genotype_file)[10:ncol(genotype_file)]){
    genotype_file[startsWith(genotype_file[[gt]], '0|0'),gt] <- 0
    genotype_file[startsWith(genotype_file[[gt]], '1|0'),gt] <- 1
    genotype_file[startsWith(genotype_file[[gt]], '0|1'),gt] <- 1
    genotype_file[startsWith(genotype_file[[gt]], '1|1'),gt] <- 2
    genotype_file[startsWith(genotype_file[[gt]], '.|.'),gt] <- NA
    genotype_file[startsWith(genotype_file[[gt]], '0/0'),gt] <- 0
    genotype_file[startsWith(genotype_file[[gt]], '1/0'),gt] <- 1
    genotype_file[startsWith(genotype_file[[gt]], '0/1'),gt] <- 1
    genotype_file[startsWith(genotype_file[[gt]], '1/1'),gt] <- 2
    genotype_file[startsWith(genotype_file[[gt]], './.'),gt] <- NA
    genotype_file[[gt]] <- as.numeric(genotype_file[[gt]])
  }
  return(genotype_file)
}

# harmonize the REF/ALT between the two genotype files
harmonized_ref_alt <- function(genotype_file1, genotype_file2){
  for(gt in colnames(genotype_file2)[10:ncol(genotype_file2)]){
    # swap homo ref and alt if necessary
    if(nrow(genotype_file2[genotype_file2$REF != genotype_file1$REF & 
                           genotype_file2$ALT != genotype_file1$ALT &
                           !is.na(genotype_file2[[gt]]) & genotype_file2[[gt]] == 0, ]) > 0){
      genotype_file2[genotype_file2$REF != genotype_file1$REF & 
                       genotype_file2$ALT != genotype_file1$ALT &
                       !is.na(genotype_file2[[gt]]) & genotype_file2[[gt]] == 0, ][[gt]] <- 2
    }
    if(nrow(genotype_file2[genotype_file2$REF != genotype_file1$REF & 
                           genotype_file2$ALT != genotype_file1$ALT &
                           !is.na(genotype_file2[[gt]]) & genotype_file2[[gt]] == 2, ]) > 0){
      genotype_file2[genotype_file2$REF != genotype_file1$REF & 
                       genotype_file2$ALT != genotype_file1$ALT &
                       !is.na(genotype_file2[[gt]]) & genotype_file2[[gt]] == 2, ][[gt]] <- 0
    }
  }
  return(genotype_file2)
}

# correlate from two genotype files and put the result in a matrix
get_gt_correlation <- function(genotype_file1, genotype_file2, verbose=T ){
  # get the participants/clusters/etc.
  geno1_cols <- colnames(genotype_file1)[10:ncol(genotype_file1)]
  geno2_cols <- colnames(genotype_file2)[10:ncol(genotype_file2)]
  # create df
  cor_matrix <- data.frame(matrix(NA, nrow=length(geno1_cols), ncol=length(geno2_cols)))
  rownames(cor_matrix) <- geno1_cols
  colnames(cor_matrix) <- geno2_cols
  # check each genotype
  for(gt1 in geno1_cols){
    for(gt2 in geno2_cols){
      # get SNPs known in both instances
      known_snps_geno1 <- genotype_file1[!is.na(genotype_file1[[gt1]]), ]$ID
      known_snps_geno2 <- genotype_file2[!is.na(genotype_file2[[gt2]]), ]$ID
      known_snps <- intersect(known_snps_geno1, known_snps_geno2)
      # grab genotypes
      geno1_gts <- genotype_file1[genotype_file1$ID %in% known_snps, ][[gt1]]
      geno2_gts <- genotype_file2[genotype_file2$ID %in% known_snps, ][[gt2]]
      # check correlation
      correlation <- cor(geno1_gts, geno2_gts, method='pearson')
      # add value in matrix
      cor_matrix[gt1, gt2] <- correlation
    }
  }
  if(verbose){
    print(cor_matrix)
  }
  return(cor_matrix)
}

# add assignments clusters, adding the best correlation matrix match to the souporcell clusters
add_assignment_to_clusters <- function(clusters, correlations, verbose=T){
  # add new columns
  clusters$assignment_ll <- NA
  clusters$correlation_ll <- NA
  # go through the clusters by grabbing the singlet assignments
  for(cluster in unique(clusters[clusters$status == 'singlet', ]$assignment)){
    # get the highest correlation
    max_cor <- max(correlations[[cluster]])
    # get the assignment
    best_cor <- rownames(correlations)[correlations[[cluster]] == max_cor]
    # add these two values to the clusters object
    clusters[clusters$assignment == cluster, ]$assignment_ll <- best_cor
    clusters[clusters$assignment == cluster, ]$correlation_ll <- max_cor
    # print assignment if necessary
    if(verbose){
      print(paste('assigning', cluster, 'to', best_cor, 'with cor of', max_cor))
    }
  }
  if(length(unique(clusters[clusters$status == 'singlet', ]$assignment)) != length(unique(clusters[clusters$status == 'singlet', ]$assignment_ll))){
    print('less correlated genotypes than clusters! sample swap?')
  }
  return(clusters)
}

# get the command line arguments
args = commandArgs(trailingOnly=TRUE)
# grab locations of genotype files
soup_vcf_loc <- args[1]
geno_vcf_loc <- args[2]
# grab of the clustering file
clus_loc <- args[3]
output_file <- args[4]

#<TEST>
# grab locations of genotype files
#soup_vcf_loc <- ''
#geno_vcf_loc <- ''
# grab of the clustering file
#clus_loc <- ''
#output_file <- ''
#</TEST>

# read the genotype files
soup_vcf <- fread(soup_vcf_loc)
geno_vcf <- fread(geno_vcf_loc)

# get the SNPs in both files
common_snps <- intersect(soup_vcf$ID, geno_vcf$ID)
# restrict to those SNPs
soup_vcf <- soup_vcf[soup_vcf$ID %in% common_snps]
geno_vcf <- geno_vcf[geno_vcf$ID %in% common_snps]

# convert to numbers
soup_vcf_conv <- convert_genotypes(soup_vcf)
geno_vcf_conv <- convert_genotypes(geno_vcf)

# harmonize soup to geno vcf alt/ref status
soup_vcf_conv_harm <- harmonized_ref_alt(geno_vcf_conv, soup_vcf_conv)

# get the correlation matrix
correlation_matrix <- get_gt_correlation(geno_vcf_conv, soup_vcf_conv_harm)

# read the clusters file
clusters <- read.table(clus_loc, header=T)

# add the new data to the clusters
clusters_correlated <- add_assignment_to_clusters(clusters, correlation_matrix)

# write the correlated clusters
write.table(clusters_correlated, output_file, quote=F, row.names=F, col.names = T, sep = '\t')
