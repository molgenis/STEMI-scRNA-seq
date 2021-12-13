#!/bin/bash
CELLRANGER_LOC='/groups/umcg-franke-scrna/tmp04/cellranger-3.1.0/cellranger'
REFDATA_LOC='/groups/umcg-franke-scrna/tmp04/refdata-cellranger-GRCh38-3.0.0/'
SEQUENCE_LANES_LOC='/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/sequence_data/'
OUTPUT_LOC='/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/alignment/b38/cellranger_output/'
OUTPUT_LOC_TMP='${TMPDIR}/releases/blokland-2020/v1/alignment/b38/cellranger_output/'
JOB_DIR='/groups/umcg-franke-scrna/tmp04/releases/blokland-2020/v1/alignment/b38/jobs/'
# these are the lanes to run through cellranger
LANES=('191125_lane1' '191125_lane2' '191126_lane1' '191126_lane2' '191209_lane1' '191209_lane2')

CORES='22'
MEMORY='190GB'
TMP_SIZE='1024gb'
RUNTIME='167:59:59'


# check each run
for lane in ${LANES[*]}
  do
    JOB_NAME='align_'${lane}
    JOB_LOC=${JOB_DIR}'/'${JOB_NAME}'_SBATCH.sh'
    JOB_OUT=${JOB_DIR}'/'${JOB_NAME}'.out'
    JOB_ERR=${JOB_DIR}'/'${JOB_NAME}'.err'

    SEQUENCE_DIR_FULL=${SEQUENCE_LANES_LOC}'/'${lane}'/'
    # temporary and permanent storage for cellranger output per lane
    OUTPUT_LOC_FULL=${OUTPUT_LOC}'/'${lane}'/'
    OUTPUT_LOC_TMP_FULL=${OUTPUT_LOC_TMP}'/'${lane}'/'
    # a library file needs to be present
    LIB_FILE=${OUTPUT_LOC_TMP_FULL}'library.csv'

    # echo the header
    echo '#!/bin/bash
#SBATCH --job-name='${JOB_NAME}'
#SBATCH --output='${JOB_OUT}'
#SBATCH --error='${JOB_ERR}'
#SBATCH --time='${RUNTIME}'
#SBATCH --cpus-per-task='${CORES}'
#SBATCH --mem='${MEMORY}'
#SBATCH --nodes=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --tmp='${TMP_SIZE}'
'> ${JOB_LOC}

    # do the prerequisites
    echo 'mkdir -p '${OUTPUT_LOC_TMP}'/'${lane}'/
cd '${OUTPUT_LOC_TMP}'
' >> ${JOB_LOC}

    # create a library file
    echo 'echo "fastqs,sample,library_type" > '${LIB_FILE} >> ${JOB_LOC}
    echo 'echo "'${SEQUENCE_DIR_FULL}','${lane}',Gene Expression" >> '${LIB_FILE} >> ${JOB_LOC}

    # build the job
    echo ${CELLRANGER_LOC}' count \
--id=cellranger_'${lane}' \
--transcriptome='${REFDATA_LOC}' \
--libraries='${LIB_FILE}' \
--localcores='${CORES}'
' >> ${JOB_LOC}

  # do the wrapping up
  echo 'cp -r '${OUTPUT_LOC_TMP_FULL}' '${OUTPUT_LOC}'
cd '${OUTPUT_LOC}'
rm -r '${OUTPUT_LOC_TMP}'
' >> ${JOB_LOC}

  done
