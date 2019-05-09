#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#                                 Running VCFANNO                                     #
#                                                                                     #
#######################################################################################

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h] [-f <vcf file>] [-p <number of processes>]
       -- script that runs vcfanno https://github.com/brentp/vcfanno/ --
       
where:
    -h    show this help text
    -s    sample name (required)
    -p    [default: 8] number of processes
"
##___________________________##
# exit script if there is no arguments
: ${1?"$usage"}

while getopts ':h:s:p:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    s) samplename=$OPTARG
       ;;
    p) nprocesses=$OPTARG
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
#Set number of processes by default
if [ -z "$nprocesses" ] ; then
  nprocesses=8
fi
##___________________________##
# Set vcf file
vcf=$(ls ${LOVELACE}Finished/${samplename}/bcbio/${samplename}/final/*/*-ensemble-annotated.vcf.gz)
##___________________________##
# Set vcfanno config file
vcfanno_config=${LOVELACE}config/grch37.vcfanno_germline.conf
##___________________________##
# Set logsample
logsample="log_${samplename}"

##################################################################
# 1) Prepare bcbio run
##################################################################
mkdir -p ${LOVELACE}Finished/${samplename}/vcfanno
cd ${LOVELACE}Finished/${samplename}/vcfanno
##################################################################
# 2) Run vcfanno
##################################################################
vcfanno -p ${nprocesses} ${vcfanno_config} ${vcf} > ${LOVELACE}Finished/${samplename}/vcfanno/${samplename}-ensemble-vcfanno.vcf

bgzip ${LOVELACE}Finished/${samplename}/vcfanno/${samplename}-ensemble-vcfanno.vcf
tabix ${LOVELACE}Finished/${samplename}/vcfanno/${samplename}-ensemble-vcfanno.vcf.gz


