##################################################################################################
## Script for  Table1 generation of STEMI patients
##################################################################################################
# Load libraries
library(Seurat)
library(ggplot2)
library(lubridate)
library(tidyr)

# Load object
clin_merged <- read.table('/Users/irene/Documents/Geneeskunde/MD-PhD/clinical_var/scRNAseq_clinical_data20201110_final_notextfields_mono_wcellnrs20201012.tsv', sep = '\t', header = T, dec = ",")

# Converting variables from factor to numneric in order to convert them
clin_merged$crp = as.numeric(as.character(clin_merged$crp))
clin_merged$peak_ck = as.numeric(as.character(clin_merged$peak_ck))
clin_merged$peak_ckmb = as.numeric(as.character(clin_merged$peak_ckmb))
clin_merged$peak_troponint = as.numeric(as.character(clin_merged$peak_troponint))
clin_merged$ischemia_time = as.numeric(as.character(clin_merged$ischemia_time))
clin_merged$bmi = as.numeric(as.character(clin_merged$bmi))
clin_merged$sys_bp = as.numeric(as.character(clin_merged$sys_bp))
clin_merged$dys_bp = as.numeric(as.character(clin_merged$dys_bp))
clin_merged$hr = as.numeric(as.character(clin_merged$hr))
clin_merged$hypertension = as.numeric(as.character(clin_merged$hypertension))
clin_merged$hypercholesterolemia = as.numeric(as.character(clin_merged$hypercholesterolemia))
clin_merged$family = as.numeric(as.character(clin_merged$family))
clin_merged$lymphocytes = as.numeric(as.character(clin_merged$lymphocytes))
clin_merged$neutrophils = as.numeric(as.character(clin_merged$neutrophils))
clin_merged$eosinophils = as.numeric(as.character(clin_merged$eosinophils))
clin_merged$basophils = as.numeric(as.character(clin_merged$basophils))
clin_merged$granulocytes_undif = as.numeric(as.character(clin_merged$granulocytes_undif))
clin_merged$monocytes = as.numeric(as.character(clin_merged$monocytes))

# Generating age variable
clin_merged$date_birth <- dmy(clin_merged$date_birth)
clin_merged$admission_date <- dmy(clin_merged$admission_date)
clin_merged$age <- interval(start= clin_merged$date_birth, end=clin_merged$admission_date)/                      
  duration(n=1, unit="years")

# Checking distributions of non-categorical variables
hist(clin_merged$age, data=clin_merged) #normal
hist(clin_merged$bmi, data=clin_merged) #normal
hist(clin_merged$dys_bp, data=clin_merged) #normal
hist(clin_merged$sys_bp, data=clin_merged) #normal
hist(clin_merged$hr, data=clin_merged) #normal
hist(clin_merged$ischemia_time, data=clin_merged) #skewed
hist(clin_merged$ck, data=clin_merged) #skewed
hist(clin_merged$ck_mb, data=clin_merged) #skewed
hist(clin_merged$troponin_t, data=clin_merged) #skewed
hist(clin_merged$nt_probnp, data=clin_merged) #skewed
hist(clin_merged$thrombocytes, data=clin_merged) #normal
hist(clin_merged$neutrophils, data=clin_merged) #normal
hist(clin_merged$lymphocytes, data=clin_merged) #skewed
hist(clin_merged$monocytes, data=clin_merged) #skewed
hist(clin_merged$eosinophils, data=clin_merged) #skewed
hist(clin_merged$basophils, data=clin_merged) #skewed
hist(clin_merged$granulocytes_undif, data=clin_merged) #skewed
hist(clin_merged$crp, data=clin_merged) #skewed
hist(clin_merged$peak_ck, data=clin_merged) #skewed
hist(clin_merged$peak_ckmb, data=clin_merged) #skewed
hist(clin_merged$peak_troponint, data=clin_merged) #skewed

# Generation of Table 1 baseline characteristics
table1(~ gender+age+bmi+hypertension+hypercholesterolemia+smoking_status+family+dys_bp+sys_bp+hr+culprit_vessel+timi_flow+timi_flow_postpci+ischemia_time+ck+ck_mb+troponin_t+nt_probnp+thrombocytes+neutrophils+lymphocytes+monocytes+eosinophils+basophils+granulocytes_undif+crp+peak_ck+peak_ckmb+peak_troponint, data=clin_merged)

# Saving file
table1(~ gender(contn\))