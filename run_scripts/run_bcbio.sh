#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#                              Running BCBIO-NEXTGEN                                  #
#                                                                                     #
#######################################################################################
#                                                                                     #
# This script runs BCBIO-NEXTGEN for a single sample. If fastq files are not defined, #
# this scripts will attempt to look for both of them (paired-end) in                  #
# ${LOVELACE}Analysis/${samplename}/fastp/                                            #
# By default, germline analysis template ${LOVELACE}config/grch37.bcbio_germline.yaml #
# will be used.                                                                       #
# Successful samples will be moved to ${LOVELACE}Finished/ and the end of the run.    #
#                                                                                     #
# REQUIRED:                                                                           #
#   - sample name [-s]                                                                #
# OPTIONAL:                                                                           #
#   - fastq1 [-f]                                                                     #
#   - fastq2 [-r]                                                                     #
#   - template [-t]                                                                   #
#                                                                                     #
####################################################################################### 

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h] [-s <sample name>] [-f <forward_1.fastq>] [-r <reverse_2.fastq>] [-t <template analysis.yaml>]
       -- script that runs bcbio https://bcbio-nextgen.readthedocs.io/en/latest/ --
       
where:
    -h    show this help text
    -s    sample name (required)
    -f    [default: samplename/fastp/*_1.fastp.fq.gz] forward fastq file
    -r    [default: samplename/fastp/*_2.fastp.fq.gz] reverse fastq file
    -t    [default: grch37.bcbio_germline.yaml] template's analysis
"
##___________________________##
# exit script if there is no arguments
: ${1?"$usage"}

while getopts ':h:s:f:r:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    t) template=$OPTARG
       ;;
    s) samplename=$OPTARG
       ;;
    f) Ffastq=$OPTARG
       ;;
    r) Rfastq=$OPTARG
       ;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

##__________ SETUP __________##
#PATHS
GENOMEDARCHIVE="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
if [ ! -d "${GENOMEDARCHIVE}" ]; then
  GENOMEDARCHIVE="/genomedarchive/"
fi
CRICK=${GENOMEDARCHIVE}Crick_storage/
LOVELACE=${GENOMEDARCHIVE}Lovelace_decoding/
MENDEL=${GENOMEDARCHIVE}Mendel_annotating/
##___________________________##
# Set foward and reverse fastq if it was not defined by user
if ([ -z "$Ffastq" ] || [ -z "$Rfastq" ]) && [ ! -d "${LOVELACE}Analysis/${samplename}/bcbio/${samplename}" ] ; then
  Ffastq=$(ls ${LOVELACE}Analysis/${samplename}/fastp/*_1.fastp.fq.gz)
  Rfastq=$(ls ${LOVELACE}Analysis/${samplename}/fastp/*_2.fastp.fq.gz)
fi
# Set germline analysis by default
template=$(ls ${LOVELACE}config/grch37.bcbio_germline.yaml)
# Set logsample
logsample="log_${samplename}"

##################################################################
# 1) Prepare bcbio run
##################################################################
# check if bcbio folder exists
# if NOT, creates it, configures BCBIO and STARTs LOG file for analysis
if [ ! -d "${LOVELACE}Analysis/${samplename}/bcbio/${samplename}" ]; then

  mkdir -p ${LOVELACE}Analysis/${samplename}/bcbio
  cd ${LOVELACE}Analysis/${samplename}/bcbio

  # Get configuration run for sample
  bcbio_nextgen.py -w template ${template} ${samplename} ${Ffastq} ${Rfastq} -n 12
  # Enter in working directory
  cd ${LOVELACE}Analysis/${samplename}/bcbio/${samplename}/work
  # Change sample name in configuration file
  sed -i "s/description:.*/description: ${samplename}/g" ../config/${samplename}.yaml

  #_________________________________
  # Start LOG file
  #_________________________________
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Starting bcbio-nextgen analysis for ${samplename}" >> ${LOVELACE}log/${logsample}
  echo "[   Software Version  ]" "bcbio-nextgen v$(bcbio_nextgen.py --version)" >> ${LOVELACE}log/${logsample}
  echo "[     File Version    ]" "${template}" >> ${LOVELACE}log/${logsample}
  echo "                       " "==> fq1: ${Ffastq}" >> ${LOVELACE}log/${logsample}
  echo "                       " "==> fq2: ${Rfastq}" >> ${LOVELACE}log/${logsample}

# if folder EXISTS, CONTINUEs LOG file for analysis
else
  #_________________________________
  # Continue LOG file
  #_________________________________
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Restarting bcbio-nextgen analysis for ${samplename}" >> ${LOVELACE}log/${logsample}
  echo "[   Software Version  ]" "bcbio-nextgen v$(bcbio_nextgen.py --version)" >> ${LOVELACE}log/${logsample}
  echo "[     File Version    ]" "${template}" >> ${LOVELACE}log/${logsample}
  echo "                       " "==> fq1: ${Ffastq}" >> ${LOVELACE}log/${logsample}
  echo "                       " "==> fq2: ${Rfastq}" >> ${LOVELACE}log/${logsample}

fi

##################################################################
# 2) Run bcbio-nextgen
##################################################################
cd ${LOVELACE}Analysis/${samplename}/bcbio/${samplename}/work
bcbio_nextgen.py ../config/${samplename}.yaml -n 12

#echo "ERROR!!"

# check if bcbio-nextgen run successfully
finished=$(grep "Timing: finished" ./log/bcbio-nextgen.log)
if [[ ! -z "$finished" ]]; then
  mv ${LOVELACE}Analysis/${samplename} ${LOVELACE}Finished/
  echo "                       " "bcbio-nextgen finished with SUCCESS for ${samplename}" >> ${LOVELACE}log/${logsample}
else
  echo "[        ERROR!!      ]" "bcbio-nextgen failed for ${samplename}" >> ${LOVELACE}log/${logsample}
  echo "                       " "This folder will not be moved, and parallel_bcbio.sh will restart the analysis as soon as possible..." >> ${LOVELACE}log/${logsample}
fi


