#!/usr/bin/env Rscript
############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_export_clinical_params.R
# Function: export clinical parameters of STEMI scRNA study for Wolf replication
############################################################################################################################

####################
# libraries        #
####################

####################
# Functions        #
####################

####################
# Main code        #
####################

# location of the sample mapping file
sample_mapping_loc <- '/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/metadata/stemi_anon_id_to_geno_id.tsv'
# and the actual clinical parameters
clinparam_loc <- '/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/clinical_parameters/scRNAseq_clinical_data_20200910_final_notextfields.tsv'
# read both tables
sample_mapping <- read.table(sample_mapping_loc, header = T, sep = ' ')
clinparam <- read.table(clinparam_loc, header = T, sep = '\t', dec = ',')
# update the record id with the TEST identifier
clinparam[['test_id']] <- paste('TEST', clinparam[['record_id']], sep = '_')
# then add the sample id
clinparam[['sample_id']] <- sample_mapping[match(clinparam[['test_id']], sample_mapping[['geno_id']]), ][['sample']]
# get the birth dates
times_births<-as.Date(clinparam[['date_birth']], format='%d/%m/%Y')
# and admission dates
times_admissions<-as.Date(clinparam[['admission_date']], format='%d/%m/%Y')
# then difference
age_diff <- difftime(times_admissions, times_births, units = 'days')
# finally as years
age_diff_years <- floor(as.vector(age_diff) / 365.25)
# add the age
clinparam[['age_at_admission']] <- age_diff_years
# now subset the clinical parameters to what we can share
clinparam_share <- clinparam[
  !is.na(clinparam[['sample_id']]), 
                             c(
                               'sample_id', 
                               'age_at_admission', 
                               'gender', 
                               'statin', 
                               'neutrophils', 
                               'lymphocytes', 
                               'non_cause', 
                               'discharge_diagnosis', 
                               'ap',
                               'mi', 
                               'pci', 
                               'cabg', 
                               'hf', 
                               'arterydisease', 
                               'tia', 
                               'stroke', 
                               'thromboembolism', 
                               'crp',
                               'crp_3h', 
                               'crp_6h', 
                               'crp_9h', 
                               'crp_12h', 
                               'crp_24h', 
                               'crp_48h'
                             )]
# aspirin is just true for all of them
clinparam_share[['aspirin']] <- 1
# save the table
clinparam_share_loc <- '/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/clinical_parameters/stemi_clinical_parameters_heinig.tsv.gz'
write.table(clinparam_share, gzfile(clinparam_share_loc), sep = '\t', quote = F, row.names = F, col.names = T)
# make a checksum
mdfiver::create_sha256_for_file(clinparam_share_loc)
