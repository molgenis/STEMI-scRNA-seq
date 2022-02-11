

#############
# libraries #
#############

library(data.table)


#############
# classes   #
#############

setClass('snp',
         slots = list(
           identifier='character',
           chrom_hg19='character',
           chrom_b38='character',
           position_hg19='numeric',
           position_b38='numeric',
           gc_score='numeric',
           gt_score='numeric',
           ref='character',
           alt='character',
           af='numeric',
           maf='numeric',
           r2='numeric',
           gt='character'
         ))


#############
# functions #
#############



report_to_vcf <- function(report, split_char='/'){
  # get the unique participants
  participants <- unique(report[['sample']])
  # get the snps
  snps <- unique(report[['SNP Name']])
  # create an empty matrix
  gt_matrix <- matrix(, ncol = length(participants), nrow = length(snps))
  # process each entry
  apply(report, 1, function(row){
    # grab the columns we care about
    snp_name <- as.vector(unlist(row[['SNP Name']]))
    sample <- as.vector(unlist(row[['Sample ID']]))
    #ilmn_strand <- row['ILMN Strand']
    #customer_strand <- row['Customer Strand']
    #allele_1_forward <- row['Allele1 - Forward']
    #allele_2_forward <- row['Allele2 - Forward']
    snp <- as.vector(unlist(row[['SNP']]))
    # get the reference and alt allele
    split_snp <- strsplit(snp, split_char)
    var_a <- split_snp[[1]][[1]]
    var_b <- split_snp[[1]][[2]]
    # get the A or the B
    allele1_ab <- as.vector(unlist(row['Allele1 - AB']))
    allele2_ab <- as.vector(unlist(row['Allele2 - AB']))
    # set up the annotation
    annotation <- '/'
    if(allele1_ab == 'A'){
      annotation <- paste('0', annotation, sep = '')
    }
    else if(allele1_ab == 'B'){
      annotation <- paste('1', annotation, sep = '')
    }
    if(allele2_ab == 'A'){
      annotation <- paste(annotation, '0', sep = '')
    }
    if(allele2_ab == 'B'){
      annotation <- paste(annotation, '1', sep = '')
    }
    # flip to be the same
    if(annotation == '1/0'){
      annotation <- '0/1'
    }
    # set in the matrix
    gt_matrix[snp, sample] <- annotation
  })
  return(report)
}


#############
# main code #
#############

# location of genotype data
genotypes_stemi_v2_report_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/genotype/Hilde_Groot_07022020_FinalReport.txt'


# read genotype data
genotypes_stemi_v2_report <- fread(genotypes_stemi_v2_report_loc, skip=9)

