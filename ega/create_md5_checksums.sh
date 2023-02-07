#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: create_md5_checksums.sh
# Function: create md5 checksums for each file
############################################################################################################################


# location of the sequence data folders
SEQUENCE_DATA_LOC='/groups/umcg-franke-scrna/prm03/releases/blokland-2020/v1/sequence_data/'
# the output for the EGA files
OUTPUT_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/ega/sequence_data/md5/'

# the lanes we want to do
LANES=('181010_lane3' '181011_lane3' '181017_lane3' '181105_lane3' \
'181106_lane3' '181121_lane3' '181122_lane3' '181213_lane4' \
'181214_lane3' '191016_lane1-ADT' '191125_lane1' '191125_lane2' \
'191126_lane1' '191126_lane2' '191209_lane1' '191209_lane2')

# check each lane
for lane in ${LANES[*]}
  do
  # create a new folder to put the md5 in
  mkdir -p ${OUTPUT_LOC}${lane}/
  # paste together the full folder
  LANE_FOLDER=${SEQUENCE_DATA_LOC}${lane}'/'
  # check each file
  for fq in ${LANE_FOLDER}*.fastq.gz
    do
    COMMAND='md5sum '${fq} > ${OUTPUT_LOC}${lane}'/'${fq##*/}'.md5'
    ${COMMAND}
  done
done
