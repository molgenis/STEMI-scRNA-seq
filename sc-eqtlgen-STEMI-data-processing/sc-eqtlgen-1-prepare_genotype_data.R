

#############
# libraries #
#############


#############
# functions #
#############

get_genotyped_locations_from_infos <- function(genotype_loc, geno_prepend, info_prepend, info_append, chromosomes=1:22, genotype_call_column='Genotyped', genotype_call_name='Genotyped', split=':', v2=T){
  # create the entire table
  full_locations_table <- NULL
  # check each of the chromosomes
  for(chrom in chromosomes){
    # paste together the location
    info_loc_chrom <- NULL
    if(v2){
      info_loc_chrom <- paste(genotype_loc, geno_prepend, chrom, '/', info_prepend, chrom, info_append, sep = '')
    }
    else{
      info_loc_chrom <- paste(genotype_loc, info_prepend, chrom, info_append, sep = '')
    }
    # read the table
    table_chrom <- read.table(info_loc_chrom, sep = '\t', header = T, stringsAsFactors = F)
    # subset to what is genotyped
    typed <- table_chrom[table_chrom[[genotype_call_column]] == genotype_call_name, 'SNP']
    # split them into parts
    typed_positions_split <- strsplit(typed, ':')
    # create an empty dataframe
    positions <- data.frame(matrix(, nrow = length(typed_positions_split), ncol = 2))
    colnames(positions) <- c('CHROM', 'POS')
    # check each entry
    for(i in 1:length(typed_positions_split)){
      # add to the dataframe
      positions[i, c('CHROM', 'POS')] <- typed_positions_split[[i]][c(1,2)]
    }
    # add to the total table
    if(is.null(full_locations_table)){
      full_locations_table <- positions
    }
    else{
      full_locations_table <- rbind(full_locations_table, positions)
    }
  }
  return(full_locations_table)
}


#############
# main code #
#############

# get the non-imputed data
# the string info for the v2 positions
v2_genotype_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/'
v2_geno_prepend <- 'chr_'
v2_info_prepend <- 'chr'
v2_info_append <- '.info.gz'
v2_positions_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v2/non_imputed_positions_hg19.tsv'
# get the infos
v2_locations <- get_genotyped_locations_from_infos(v2_genotype_loc, v2_geno_prepend, v2_info_prepend, v2_info_append)
# write the result
write.table(v2_locations, v2_positions_output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
# subset to exising data using position file
#bcftools view -O z -o ./stemi_v2_unimputed.vcf.gz -R ../non_imputed_positions_hg19.tsv /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20190411/genotypes/allTestChroms.vcf.gz
# annotate with RSids
#bcftools annotate -O z -o ./stemi_v2_unimputed.rsid.vcf.gz -a /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz -c ID stemi_v2_unimputed.vcf.gz
# harmonize the ref/alt positions
#singularity exec --bind /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v2/genotype/unimputed/,/groups/umcg-franke-scrna/tmp01/external_datasets/,/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/ /groups/umcg-wijmenga/tmp01/users/umcg-roelen/singularity/single_cell_container/singlecell_container.simg bcftools +fixref ./stemi_v2_unimputed.rsid.vcf.gz -O z -o ./stemi_v2_unimputed.rsid.reffixed.vcf.gz -- -d -f /groups/umcg-franke-scrna/tmp01/external_datasets/refdata-cellranger-hg19-3.0.0/fasta/genome.fa -i /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz

# the string info for the v3 positions
v3_genotype_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20200207/imputed/'
v3_geno_prepend <- NA
v3_info_prepend <- 'chr'
v3_info_append <- '.info.gz'
v3_positions_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/non_imputed_positions_hg19.tsv'
# get the infos
v3_locations <- get_genotyped_locations_from_infos(v3_genotype_loc, v3_geno_prepend, v3_info_prepend, v3_info_append, v2 = F)
# write the result
write.table(v3_locations, v3_positions_output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
# subset the existing genotype data using these position files
# bcftools view -O z -o ./stemi_v3_unimputed.vcf.gz -R ../non_imputed_positions_hg19.tsv /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20200207/imputed/all_chr.vcf.gz
v3_genotype_pt88_loc <- '/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20200617/TEST_88/'
v3_positions_pt_88_output_loc <- '/groups/umcg-franke-scrna/tmp04/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/non_imputed_positions_hg19_pt88.tsv'
v3_locations_pt88 <- get_genotyped_locations_from_infos(v3_genotype_pt88_loc, v3_geno_prepend, v3_info_prepend, v3_info_append, v2 = F)
write.table(v3_locations_pt88, v3_positions_pt_88_output_loc, sep = '\t', row.names = F, col.names = T, quote = F)
# subset the existing genotype data using these position files
#bcftools view -O z -o ./stemi_v3_pt88_unimputed.vcf.gz -R ../non_imputed_positions_hg19_pt88.tsv /groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/20200617/TEST_88/chrAll.TEST_88.vcf.gz
# annotate with RSids
#bcftools annotate -O z -o ./stemi_v3_unimputed.rsid.vcf.gz -a /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz -c ID stemi_v3_unimputed.vcf.gz
#bcftools annotate -O z -o ./stemi_v3_pt88_unimputed.rsid.vcf.gz -a /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz -c ID stemi_v3_pt88_unimputed.vcf.gz
# harmonize the ref/alt positions
#singularity exec --bind /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/genotype/unimputed/,/groups/umcg-franke-scrna/tmp01/external_datasets/,/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/ /groups/umcg-wijmenga/tmp01/users/umcg-roelen/singularity/single_cell_container/singlecell_container.simg bcftools +fixref ./stemi_v3_unimputed.rsid.vcf.gz -O z -o ./stemi_v3_unimputed.rsid.reffixed.vcf.gz -- -d -f /groups/umcg-franke-scrna/tmp01/external_datasets/refdata-cellranger-hg19-3.0.0/fasta/genome.fa -i /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz
#singularity exec --bind /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/genotype/unimputed/,/groups/umcg-franke-scrna/tmp01/external_datasets/,/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/ /groups/umcg-wijmenga/tmp01/users/umcg-roelen/singularity/single_cell_container/singlecell_container.simg bcftools +fixref ./stemi_v3_pt88_unimputed.rsid.vcf.gz -O z -o ./stemi_v3_pt88_unimputed.rsid.reffixed.vcf.gz -- -d -f /groups/umcg-franke-scrna/tmp01/external_datasets/refdata-cellranger-hg19-3.0.0/fasta/genome.fa -i /groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/dbSNP/human_9606_b151_GRCh37p13/All_20180423.vcf.gz

# merge the v2, v3 and pt88 vcf files
#singularity exec --bind /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/genotype/unimputed/,/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v2/genotype/unimputed/,/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/ /groups/umcg-wijmenga/tmp01/users/umcg-roelen/singularity/single_cell_container/singlecell_container.simg bcftools merge -O z -o /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/stemi_all_unimputed.vcf.gz /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v2/genotype/unimputed/stemi_v2_unimputed.rsid.reffixed.vcf.gz /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/genotype/unimputed/stemi_v3_unimputed.rsid.reffixed.vcf.gz /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_v3/genotype/unimputed/stemi_v3_pt88_unimputed.rsid.reffixed.vcf.gz
#singularity exec --bind /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/ /groups/umcg-wijmenga/tmp01/users/umcg-roelen/singularity/single_cell_container/singlecell_container.simg plink2 --vcf /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/stemi_all_unimputed.vcf.gz --make-pgen --out /groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/stemi_all_unimputed

# create the metadata files
# metadata file location
metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/metadata/stemi_age_gender_match_wtestid.tsv'
# read the metadata file
metadata <- read.table(metadata_loc, sep = '\t', header = T)

# create the required psam file
all_original_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/stemi_all_unimputed.psam.original'
all_new_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_blokland2020/stemi_all/genotype/unimputed/stemi_all_unimputed.psam'
all_psam <- read.table(all_original_psam_loc, header = F, sep = '\t')
colnames(all_psam) <- c('#FID', 'SEX')
# we remove the sex column, as it is empty right now
all_psam[['SEX']] <- NULL
# add the other ID column
all_psam[['IID']] <- all_psam[['#FID']]
# needs to be numeric?
all_psam[['PAT']] <- 0
all_psam[['MAT']] <- 0
# grab the sex from the metadata
all_psam_sex <- metadata[match(all_psam[['IID']], metadata[['ID']]), 'gender']
# create the empty sex column in the psam
all_psam[['SEX']] <- NA
# we need to change the coding from M/F to 1/2 in the psam
all_psam[!is.na(all_psam_sex) & all_psam_sex == 'M', 'SEX'] <- 1
all_psam[!is.na(all_psam_sex) & all_psam_sex == 'F', 'SEX'] <- 2
all_psam[is.na(all_psam_sex), 'SEX' ] <- 0
# we don't know most of these
all_psam[['Provided_Ancestry']] <- 'EUR'
all_psam[['genotyping_platform']] <- 'cytoSNP'
all_psam[['array_available']] <- 'N'
all_psam[['wgs_available']] <- 'N'
all_psam[['wes_available']] <- 'Y'
all_psam[['age']] <- metadata[match(all_psam[['IID']], metadata[['ID']]), 'age']
all_psam[['age_range']] <- apply(all_psam, 1, function(x){
  floor(as.numeric(x['age'])/10)*10
})
all_psam[['Study']] <- 'stemi_all'
all_psam[['smoking_status']] <- NA
all_psam[['hormonal_contraception_use_currently']] <- NA
all_psam[['menopause']] <- NA
all_psam[['pregnancy_status']] <- NA
write.table(all_psam, all_new_psam_loc, sep = '\t', row.names = F, col.names = T, quote = F)



