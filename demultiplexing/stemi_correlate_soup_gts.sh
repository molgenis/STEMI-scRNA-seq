#!/usr/bin/env bash

#directory and file listings
JOB_FILE='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/correlate_clusters/jobs/correlate.out'
COR_SCRIPT_LOC='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/correlate_clusters/scripts/stemi_correlate_soup_gts.R'
GENOTYPES_PREPEND='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/'
GENOTYPES_APPEND='_mmaf005_exons_sorted.vcf'
SOUP_GENO_PREPEND='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/souporcell_output/'
SOUP_GENO_APPEND='cluster_genotypes.vcf'
SOUP_CLUSTERS_PREPEND='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/souporcell_output/'
SOUP_CLUSTERS_APPEND='clusters.tsv'
COR_OUTPUT_PREPEND='/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/correlate_clusters/correlated_output/'
COR_OUTPUT_APPEND='_correlated.tsv'

# lanes in the project
TEST_LANES=("181010_lane3" "181011_lane3" "181017_lane3" "181105_lane3" "181106_lane3" "181121_lane3" "181122_lane3" "181213_lane4" "181214_lane3" "191125_lane1" "191125_lane2" "191126_lane1" "191126_lane2" "191209_lane1" "191209_lane2")
#TEST_LANES=("181010_lane3" "181011_lane3" "181017_lane3" "181105_lane3" "181106_lane3" "181121_lane3" "181122_lane3" "181213_lane4" "181214_lane3" "191125_lane1" "191125_lane2" "191126_lane1" "191126_lane2" "191209_lane1")

# we need R
ml R

# check each lane
for dir in ${TEST_LANES[*]} ; do

# write the lane we're at
echo ${dir} >> ${JOB_FILE}
# create full paths
genotype_loc=${GENOTYPES_PREPEND}${dir}${GENOTYPES_APPEND}
soup_geno_loc=${SOUP_GENO_PREPEND}${dir}/${SOUP_GENO_APPEND}
soup_clus_loc=${SOUP_CLUSTERS_PREPEND}${dir}/${SOUP_CLUSTERS_APPEND}
cor_output=${COR_OUTPUT_PREPEND}${dir}${COR_OUTPUT_APPEND}

# do correlation and write any logging to job file
Rscript ${COR_SCRIPT_LOC} \
${soup_geno_loc} \
${genotype_loc} \
${soup_clus_loc} \
${cor_output} \
>> ${JOB_FILE}

done
