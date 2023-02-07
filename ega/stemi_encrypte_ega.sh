#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_encrypt_ega.sh
# Function: encrypt the files for upload to EGA
############################################################################################################################


# location of the sequence data folders
SEQUENCE_DATA_LOC='/groups/umcg-franke-scrna/prm03/releases/blokland-2020/v1/sequence_data/'
# the output for the EGA files
OUTPUT_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/ega/sequence_data/'

# the lanes we want to do
LANES=('181010_lane3' '181011_lane3' '181017_lane3' '181105_lane3' \
'181106_lane3' '181121_lane3' '181122_lane3' '181213_lane4' \
'181214_lane3' '191016_lane1-ADT' '191125_lane1' '191125_lane2' \
'191126_lane1' '191126_lane2' '191209_lane1' '191209_lane2')

# check each lane
for lane in ${LANES[*]}
  do

  # paste together the full folder
  LANE_FOLDER=${SEQUENCE_DATA_LOC}${lane}'/'
  # check each file
  for fq in ${LANE_FOLDER}*.fastq.gz
    do
    # create the command
    ENCRYPT_COMMAND="java -jar EGA-Cryptor-2.0.0/ega-cryptor-2.0.0.jar -i "${fq}" -o "${OUTPUT_LOC}
    # and do the encryption
    ${ENCRYPT_COMMAND}
  done
done
