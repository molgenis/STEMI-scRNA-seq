#!/usr/bin/env bash

#directory and file listings
LANES_DIR="/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/cellranger_output/"
LANE_READGROUP_APPEND="outs/possorted_genome_bam.bam"
LANE_BARCODE_APPEND="outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
INDIVIDUAL_GENOTYPES="/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/genotype/genotype_per_lane_mmaf005_exons/"
GENOTYPE_APPEND="_mmaf005_exons_sorted.vcf"
OUTPUT_DIR="/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/souporcell_output/"
JOB_DIR="/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/demultiplexing/souporcell/jobs/"
SOUPOR_IMAGE="/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/image/souporcell.sif"
GENOME_LOC="/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/external_software/refdata-cellranger-hg19-3.0.0/fasta/genome.fa"

TEST_LANES=("181010_lane3" "181011_lane3" "181017_lane3" "181105_lane3" "181106_lane3" "181121_lane3" "181122_lane3" "181213_lane4" "181214_lane3" "191125_lane1" "191125_lane2" "191126_lane1" "191126_lane2" "191209_lane1" "191209_lane2")

for dir in ${TEST_LANES[*]} ; do

    lane_id=${dir}
    output_folder=${OUTPUT_DIR}/${lane_id}/
    mkdir -p ${output_folder}
    output_job=${JOB_DIR}/soup_${lane_id}_stemi_mmaf005_SBATCH.sh

    echo -e "#!/usr/bin/env bash

#SBATCH --job-name=soup_${lane_id}_mmaf005
#SBATCH --output=soup_${lane_id}_mmaf005.out
#SBATCH --error=soup_${lane_id}_mmaf005.err
#SBATCH --time=23:59:59
#SBATCH --cpus-per-task=8
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

        set -e

        cd ${output_folder}

        export SINGULARITY_BINDPATH=\"/groups/umcg-wijmenga/tmp01/projects/1M_cells_scRNAseq/ongoing/Cardiology/,/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/demultiplexing/souporcell/image/,/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/external_software/\"

        singularity exec ${SOUPOR_IMAGE} souporcell_pipeline.py \
-i ${LANES_DIR}/${lane_id}/outs/possorted_genome_bam.bam \
-b ${LANES_DIR}/${lane_id}/outs/filtered_feature_bc_matrix/barcodes.tsv \
-f ${GENOME_LOC} \
-t 8 \
-o ${output_folder}/ \
-k 8 \
--known_genotypes ${INDIVIDUAL_GENOTYPES}/${lane_id}${GENOTYPE_APPEND} \
--skip_remap True


        " >> ${output_job}
done
