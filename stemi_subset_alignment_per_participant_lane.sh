#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_subset_alignment_per_participant_lane.sh
# Function: create md5 checksums for each file
############################################################################################################################

# location of the binary to execute the subset
SUBSET_BAM_BIN_LOC='/groups/umcg-franke-scrna/tmp01/software/subset-bam/subset-bam_linux'
# location of the text files with subsets
SAMPLE_TXTS_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/variant_calling/barcodes_per_sample/'
# the filename append of the barcode lists
SAMPLE_TXTS_APPEND='_barcodes.txt'
# location of the output bams
SUBSET_BAMS_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/variant_calling/aligned_reads_per_lane/'
# location of the source bams
CELLRANGER_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/alignment/b38/cellranger_output/'
# path of bam file in the cellranger folders
IN_BAM_PATH='outs/possorted_genome_bam.bam'

# the lanes we want to do
LANES=('181010_lane3' '181011_lane3' '181017_lane3' '181105_lane3' \
'181106_lane3' '181121_lane3' '181122_lane3' '181213_lane4' \
'181214_lane3' '191016_lane1-ADT' '191125_lane1' '191125_lane2' \
'191126_lane1' '191126_lane2' '191209_lane1' '191209_lane2')

# check each lane
for lane in ${LANES[*]}
  do
    # get the full path to the input bam
    IN_BAM_FULL_LOC=${CELLRANGER_LOC}${lane}'/'${IN_BAM_PATH}
    # paste together the location of the sequence file
    BARCODE_FILE_PATTERN=${SAMPLE_TXTS_LOC}${lane}'/*'${SAMPLE_TXTS_APPEND}
    # match
    BARCODE_FILES_MATCH=( ${BARCODE_FILE_PATTERN} )
    # now check each file
    for barcode_file in ${BARCODE_FILES_MATCH[*]}
        do
        # extract the sample name
        BASE_FILE=$(basename -- ${barcode_file})
        # get the samplename
        SAMPLE=${BASE_FILE/${SAMPLE_TXTS_APPEND}/""}
        # make the directory we will write to
        mkdir -p ${SUBSET_BAMS_LOC}${lane}'/'
        # get the output bam full location
        OUT_BAM_FULL_LOC=${SUBSET_BAMS_LOC}${lane}'/'${SAMPLE}'.bam'
        # build the command
        SUBSET_COMMAND=${SUBSET_BAM_BIN_LOC}' --bam '${IN_BAM_FULL_LOC}' --cell-barcodes '${barcode_file}' --out-bam '${OUT_BAM_FULL_LOC}
        # do the actual command
        ${SUBSET_COMMAND}
    done
done
