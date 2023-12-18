#!/bin/bash

###################################################################
#Script Name	  : stemi_limix_input_files.sh
#Description	  : create Snakemake input files
#Args           : 
#Author       	: Roy Oelen
# example        : ./stemi_limix_input_files.sh \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/genotype/stemi_all_nc2022_biallelic \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v3/smf.txt \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v3/ \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/output/cell_type_lowerres/v3/unconfined \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/snakemake/cell_type_lowerres/v3/unconfined
#
# ./stemi_limix_input_files.sh \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/genotype/stemi_all_nc2022_aragam2022 \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v3/smf.txt \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/input/cell_type_lowerres/v3/ \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/output/cell_type_lowerres/v3/aragam2022 \
# /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/snakemake/cell_type_lowerres/v3/aragam2022
#
###################################################################

# location of the templates
SNAKEMAKE_TEMPLATE_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/templates/Snakemake.smk.py.template'
JSON_TEMPLATE_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/templates/cluster.json.template'
# the covariate file name
COVAR_FILENAME='covariates_status'

# the command parameters
GENO_LOC=$1
SAMPLE_MAPPING_LOC=$2
EXP_LOC=$3
OUTPUT_LOC=$4
SNAKEMAKE_LOC=$5


# create the regex for the seurat objects
REGEX_EXP_FILES=${EXP_LOC}'/*.Exp.txt'

# check each expression file
for exp_file in ${REGEX_EXP_FILES}
    do
    # fetch the basename
    file_basename=$(basename "${exp_file}" .Exp.txt)
    # get the expression full filename
    phenotype_loc=${EXP_LOC}'/'${file_basename}'.Exp.txt'
    # get the covariate full filename
    covariate_loc=${EXP_LOC}'/'${file_basename}'.'${COVAR_FILENAME}'.txt'
    # get the output directory
    output_dir=${OUTPUT_LOC}'/'${file_basename}'/'
    # create the directory
    mkdir -p ${output_dir}
    # get the snakemake directory
    snakemake_dir=${SNAKEMAKE_LOC}'/'${file_basename}
    # create the directory
    mkdir -p ${snakemake_dir}
    # get the snakemake loc
    snakemake_loc=${snakemake_dir}'/Snakemake.py.smk'
    # now start replacing some  
    sed 's,\[\[genotype_template_placeholder\]\],'"${GENO_LOC}"',' ${SNAKEMAKE_TEMPLATE_LOC} > ${snakemake_loc}
    sed -i 's,\[\[covariate_template_placeholder\]\],'"${covariate_loc}"',' ${snakemake_loc}
    sed -i 's,\[\[phenotype_template_placeholder\]\],'"${phenotype_loc}"',' ${snakemake_loc}
    sed -i 's,\[\[sample_mapping_template_placeholder\]\],'"${SAMPLE_MAPPING_LOC}"',' ${snakemake_loc}
    sed -i 's,\[\[output_template_placeholder\]\],'"${output_dir}"',' ${snakemake_loc}
    # same for the cluster.json
    json_loc=${snakemake_dir}'/'${file_basename}'.cluster.json'
    sed 's,\[\[cell_type\]\],'"${file_basename}"',' ${JSON_TEMPLATE_LOC} > ${json_loc}
    sed -i 's,\[\[slurm_out_dir\]\],'"${output_dir}"',' ${json_loc}
    sed -i 's,\[\[slurm_err_dir\]\],'"${output_dir}"',' ${json_loc}
done

# the snakemake pipeline is started like this:
# conda activate limix_qtl
# snakemake \
# -j 99 \
# -d /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/snakemake/cell_type_lowerres/v2/unconfined/monocyte/ \
# -s /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/snakemake/cell_type_lowerres/v2/unconfined/monocyte/Snakemake.py.smk \
# --cluster-config /groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/eqtl/sc-eqtlgen/snakemake/cell_type_lowerres/v2/unconfined/monocyte/monocyte.cluster.json \
# --latency-wait 120 \
# --cluster "sbatch -n {cluster.n} -t {cluster.time} -o {cluster.output} -e {cluster.error} --mem {cluster.memory}"
