#!/usr/bin/env bash
geno_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/stemi_all_merged_nomissing_numeric.tsv'
snppos_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/snp_pos.tsv'
genepos_loc='/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/eQTL_mapping/metadata/gene_positions.tsv'

features_loc=$1
result_loc=$2
metadata_loc=$3
job_loc=$4
job_name=$5

echo "#!/usr/bin/env bash
#SBATCH --job-name=${job_name}
#SBATCH --output=${job_name}.out
#SBATCH --error=${job_name}.err
#SBATCH --time=05:59:59
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

ml R
Rscript \
${features_loc} \
${geno_loc} \
${snppos_loc} \
${genepos_loc} \
${metadata_loc}" > ${job_loc}
fi
