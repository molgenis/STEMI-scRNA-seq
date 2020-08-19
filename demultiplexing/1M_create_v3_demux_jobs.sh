#!/usr/bin/env bash

#directory and file listings
LANES_DIR="/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/cellranger_output/"
LANE_READGROUP_APPEND="outs/possorted_genome_bam.bam"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
INDIVIDUAL_GENOTYPES="/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/"
GENOTYPE_APPEND="_mmaf005_exons_sorted.vcf"
OUTPUT_DIR="/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/demuxlet/demuxlet_output/"
JOB_DIR="/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/demuxlet/jobs/"
DEMUXLET_DIR="/groups/umcg-wijmenga/tmp04/umcg-hbrugge/apps/demuxlet/"

#parameters used
TAG_GROUP="CB"
TAG_UMI="UB"
FIELD="GT"

TEST_LANES=("191125_lane1" "191125_lane2" "191126_lane1" "191126_lane2" "191209_lane1" "191209_lane2")

for dir in ${TEST_LANES[*]} ; do

    lane_id=${dir}
    output_job=${JOB_DIR}/demux_${lane_id}_stemi_mmaf005_SBATCH.sh

    echo -e "#!/usr/bin/env bash

#SBATCH --job-name=demux_${lane_id}_mmaf005
#SBATCH --output=demux_${lane_id}_mmaf005.out
#SBATCH --error=demux_${lane_id}_mmaf005.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=1
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

        set -e

        ${DEMUXLET_DIR}demuxlet \\
                --sam ${LANES_DIR}${lane_id}/${LANE_READGROUP_APPEND} \\
                --tag-group ${TAG_GROUP} \\
                --tag-UMI ${TAG_UMI} \\
                --field ${FIELD} \\
                --vcf ${INDIVIDUAL_GENOTYPES}${lane_id}${GENOTYPE_APPEND} \\
                --out ${OUTPUT_DIR}${lane_id}_mmaf002 \\
                --group-list ${LANES_DIR}${lane_id}/${LANE_BARCODE_APPEND}

        " >> ${output_job}
done
