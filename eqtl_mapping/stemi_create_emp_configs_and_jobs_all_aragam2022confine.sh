#!/usr/bin/env bash
###################################################################
#Script Name	  : stemi_create_emp_configs_and_jobs_all_aragam2022confine.sh
#Description	  : create SBATCH jobs scripts to do EMP QTL mapping
#Args           :
#Author       	: Roy Oelen
#example        : ./stemi_create_emp_configs_and_jobs_all_aragam2022confine.sh \
#                 /groups/umcg-weersma/tmp01/projects/lpmc_v2/ongoing/qtl_mapping/expression_files/emp/elmentaite_adult_martin_immune_20220824/sct/count/AI/ \
#                 /groups/umcg-weersma/tmp01/projects/lpmc_v2/ongoing/qtl_mapping/emp/configs/elmentaite_adult_martin_immune_20220824/sct/count/AI/ \
#                 /groups/umcg-weersma/tmp01/projects/lpmc_v2/ongoing/qtl_mapping/emp/results/elmentaite_adult_martin_immune_20220824/sct/count/AI/ \
#                 /groups/umcg-weersma/tmp01/projects/lpmc_v2/ongoing/qtl_mapping/emp/jobs/elmentaite_adult_martin_immune_20220824/sct/count/AI/
###################################################################


FEATURES_LOC=$1
CONFIGS_LOC=$2
RESULTS_LOC=$3
JOBS_LOC=$4

CONFINEMENT='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/qtl_mapping/confinements/aragam_2022_primary_gws_gwas_snps.txt'
ANNOTATION_LOC='/groups/umcg-franke-scrna/tmp01/external_datasets/emp_annotations/singleCell-annotation-stripped.tsv'
EMP_LOC='/groups/umcg-franke-scrna/tmp01/software/eqtl-mapping-pipeline-1.4.0-SNAPSHOT/eqtl-mapping-pipeline.jar'
GENOTYPE_LOC='/groups/umcg-franke-scrna/tmp01/releases/blokland-2020/v1/genotype/stemi_all_merged_trityper/'
CISTRANS='cis'

cd ${FEATURES_LOC}
for file in *; do mv "$file" `echo $file | tr ' ' '_'` ; done

# loop through files
files=${FEATURES_LOC}"/*"
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
    #mv $f "${FEATURES_LOC}/${posix_basename}.tsv"
    # echo ${posix_basename}
    # create location for config file
    config_location=${CONFIGS_LOC}"/"${posix_basename}"_emp_config.xml"
    # create the folder if non_existant
    mkdir -p ${CONFIGS_LOC}
    echo ${config_location}
    # create location for the result files
    result_location=${RESULTS_LOC}"/"${posix_basename}"/"
    # create the folder if non_existant
    mkdir -p ${result_location}
    echo ${result_location}
    # create location for job file
    job_location=${JOBS_LOC}"/"${posix_basename}"_emp_job_SBATCH.sh"
    # create the folder if non_existant
    mkdir -p ${JOBS_LOC}
    echo ${job_location}

  # create the config file
echo "<?xml version=\"1.0\" encoding=\"utf-8\" standalone=\"no\"?>

<settings>

  <defaults>

    <qc>
        <snpqccallratethreshold>0.95</snpqccallratethreshold>
        <snpqchwethreshold>0.0001</snpqchwethreshold>
        <snpqcmafthreshold>0.05</snpqcmafthreshold>
    </qc>

    <analysis>
        <analysistype>${CISTRANS}</analysistype>
        <cisanalysisprobedistance>100000</cisanalysisprobedistance>
        <correlationtype>nonparametric</correlationtype>
        <threads>10</threads>
        <createdotplot>false</createdotplot>
        <createqqplot>false</createqqplot>
    </analysis>

    <multipletesting>
        <type>fdr</type>
        <threshold>0.05</threshold>
         <fdrtype>probe-level</fdrtype>
        <permutations>10</permutations>
    </multipletesting>

    <output>
        <outputdirectory>${result_location}</outputdirectory>
        <outputplotthreshold>0</outputplotthreshold>
        <outputplotdirectory>${result_location}</outputplotdirectory>
        <maxnreqtlresults>50000000</maxnreqtlresults>
        <generatesnpsummarystatistics>false</generatesnpsummarystatistics>
        <generateeqtlpvaluetable>false</generateeqtlpvaluetable>
        <binaryoutput>false</binaryoutput>
        <textoutput>true</textoutput>
    </output>

    <confine>
        <snp>${CONFINEMENT}</snp>
<snpProbe></snpProbe>
        <confineSNPsToSNPsPresentInAllDatasets>true</confineSNPsToSNPsPresentInAllDatasets>
        <confineSNPsSelectSNPInStrongestLD>false</confineSNPsSelectSNPInStrongestLD>
        <confineProbesThatMapToKnownChromosome>true</confineProbesThatMapToKnownChromosome>
    </confine>

  </defaults>

  <datasets>
    <dataset>
        <name>${posix_basename}_all</name>
		<location>${GENOTYPE_LOC}</location>
		<expressiondata>${f}</expressiondata>
		<probeannotation>${ANNOTATION_LOC}</probeannotation>
		<quantilenormalize>false</quantilenormalize>
		<logtranform>false</logtranform>
    </dataset>
  </datasets>

</settings>" > ${config_location}

# create the job file
echo "#!/usr/bin/env bash
#SBATCH --job-name=map_${posix_basename}
#SBATCH --output=${JOBS_LOC}/map_${posix_basename}.out
#SBATCH --error=${JOBS_LOC}/map_${posix_basename}.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=10
#SBATCH --mem=64gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

set -e
ml Java

java -jar -Xmx40g -Xms20g -XX:StringTableSize=10000019 -XX:MaxPermSize=512m \
${EMP_LOC} \
--mode metaqtl \
--settings ${config_location}" > ${job_location}
  fi
done
