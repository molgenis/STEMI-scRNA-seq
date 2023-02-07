#!/bin/bash

############################################################################################################################
# Authors: Roy Oelen
# Name: stemi_create_pairing_file.sh
# Function: create md5 checksums for each file
############################################################################################################################


# location of the sequence data folders
SEQUENCE_DATA_LOC='/groups/umcg-franke-scrna/prm03/releases/blokland-2020/v1/sequence_data/'
# location of the md5s that are unencrypted
MD5LOC_UNENCRYPTED='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/ega/sequence_data/md5/'
# location of the encrypted md5s
MD5LOC_ENCRYPTED='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/ega/sequence_data/'
# location of the sample mapping file
PAIRING_FILE_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/ega/sequence_data/pairing_file_lanes.csv'

# the lanes we want to do
LANES=('181010_lane3' '181011_lane3' '181017_lane3' '181105_lane3' \
'181106_lane3' '181121_lane3' '181122_lane3' '181213_lane4' \
'181214_lane3' '191016_lane1-ADT' '191125_lane1' '191125_lane2' \
'191126_lane1' '191126_lane2' '191209_lane1' '191209_lane2')

# write the header
echo 'Sample Alias,First Fastq File,First Checksum,First Unencrypted checksum,Second Fastq File,Second Checksum,Second Unencrypted checksum' > ${PAIRING_FILE_LOC}
# check each lane
for lane in ${LANES[*]}
  do
  # paste together the location of the sequence file
  SEQUENCE_R1_PATTERN=${SEQUENCE_DATA_LOC}${lane}'/*R1*.fastq.gz'
  # match
  R1_FILES_MATCH=( ${SEQUENCE_R1_PATTERN} )
  # now check each file
  for R1 in ${R1_FILES_MATCH[*]}
    do
    # get the corresponding R2
    R2=${R1/R1/R2}
    # get the basenames of the R1 and R2
    R1_NOPATH=${R1##*/}
    R2_NOPATH=${R2##*/}
    # get the location of the corresponding md5s
    R1_MD5_UNENCRYPTED_LOC=${MD5LOC_UNENCRYPTED}${lane}/${R1_NOPATH}'.md5'
    R2_MD5_UNENCRYPTED_LOC=${MD5LOC_UNENCRYPTED}${lane}/${R2_NOPATH}'.md5'
    R1_MD5_ENCRYPTED_LOC=${MD5LOC_ENCRYPTED}${R1_NOPATH}'.gpg.md5'
    R2_MD5_ENCRYPTED_LOC=${MD5LOC_ENCRYPTED}${R2_NOPATH}'.gpg.md5'
    # get the md5s from the files
    R1_UNENCRYPTED_MD5=$(sed -n '1p' ${R1_MD5_UNENCRYPTED_LOC} | awk '{print $1}' | tr -d '\n') # this file is a bit messy, we need the first row and column value [0,0] and then need to remove the newline
    R2_UNENCRYPTED_MD5=$(sed -n '1p' ${R2_MD5_UNENCRYPTED_LOC} | awk '{print $1}' | tr -d '\n')
    R1_ENCRYPTED_MD5=$(cat ${R1_MD5_ENCRYPTED_LOC})
    R2_ENCRYPTED_MD5=$(cat ${R1_MD5_ENCRYPTED_LOC})
    # now strip extention and read from the filenames so we get just the sample names, then remove the newline
    SAMPLE=$( echo ${R1_NOPATH%%_R1_*.*} | tr -d '\n' )
    # start writing to file
    # Sample Alias,First Fastq File,First Checksum,First Unencrypted checksum,Second Fastq File,Second Checksum,Second Unencrypted checksum
    LINE_TO_WRITE=${SAMPLE}','${R1_NOPATH}','${R1_ENCRYPTED_MD5}','${R1_UNENCRYPTED_MD5}','${R2_NOPATH}','${R2_ENCRYPTED_MD5}','${R2_UNENCRYPTED_MD5}
    echo ${LINE_TO_WRITE} >> ${PAIRING_FILE_LOC}
  done
done
