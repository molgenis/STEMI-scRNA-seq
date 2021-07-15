#!/usr/bin/env bash
geno_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
snppos_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
genepos_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'

#./stemi_create_matrixqtl_jobs.sh \
#/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/features/stemi_all_lowerres_20210301_metaqtl/Baseline/ \
#/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/results/MatrixeQTL/stemi_all_lowerres_20210301/Baseline/ \
#/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/stemi_all_lowerres_20210301/Baseline/ \
#/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/jobs/MatrixEQTL/stemi_all_lowerres_20210301_unconfined/Baseline/

features_loc=$1
results_loc=$2
metadata_loc=$3
jobs_loc=$4

cd ${features_loc}
for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done

# loop through files
files=${features_loc}"/*"
for f in ${files}
do
  # only create files for features
  if [[ ${f} =~ .*\.(tsv) ]]
  then
    # ${f} is the features file full path, we also need just the name
    # echo "${f}"
    fbasename=$(basename "${f}" .tsv)
    # echo "${fbasename}"
    # make basename posix compliant (i.e. no spaces etc.)
    posix_basename=${fbasename//[^a-zA-Z0-9]/_}
    #posix_basename="echo ${fbasename} | tr ' ' '_'"
    #mv $f "${features_loc}/${posix_basename}.tsv"
    # echo ${posix_basename}
    # create location for the result files
    result_location=${results_loc}"/"${posix_basename}"_matrixeqtlout.tsv"
    # create the folder if non_existant
    mkdir -p ${result_location}
    echo ${result_location}
    # create location for job file
    job_location=${jobs_loc}"/"${posix_basename}"_matrixeqtl_job_SBATCH.sh"
    # create the folder if non_existant
    mkdir -p ${jobs_loc}
    metadata_file="${posix_basename/_expression/_metadata.chem.tsv}"
    metadata_location=${metadata_loc}${metadata_file}

echo "#!/usr/bin/env bash
#SBATCH --job-name=map_${posix_basename}
#SBATCH --output=${jobs_loc}/map_${posix_basename}.out
#SBATCH --error=${jobs_loc}/map_${posix_basename}.err
#SBATCH --time=05:59:59
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml R
Rscript \
${f} \
${geno_loc} \
${snppos_loc} \
${genepos_loc} \
${result_location} \
${metadata_location}" > ${job_location}
  fi
done
