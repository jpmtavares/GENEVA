#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#                      Parallelization of BCBIO-NEXTGEN                               #
#                                                                                     #
#######################################################################################
#                                                                                     #
# This script parallelizes BCBIO-NEXTGEN for all samples in ${LOVELACE}Analysis/,     #
# using run_bcbio.sh                                                                  # 
# By default, this script will run 4 samples simultaneously by order them from older  #
# to newer. All samples that contain a bcbio folder, will be ignore in the first run, #
# because it means that they are already being runned elsewhere. They will however be #
# listed in a second run, because it means that they were runned but failed.          #
#                                                                                     #
# REQUIRED:                                                                           #
#   - run_bcbio.sh                                                                    #
# OPTIONAL:                                                                           #
#   - number of samples to run simultaneously [-n]                                    #
#                                                                                     #
#######################################################################################

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h] [-n <number of samples>]
       -- script that parallelizes bcbio https://bcbio-nextgen.readthedocs.io/en/latest/
          through all the samples present in {LOVELACE}Analysis/ --
       
where:
    -h    show this help text
    -n    number of samples that will run simultaneously
"

##__________ SETUP __________##
#PATHS
GENOMEDARCHIVE="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
if [ ! -d "${GENOMEDARCHIVE}" ]; then
  GENOMEDARCHIVE="/genomedarchive/"
fi
CRICK=${GENOMEDARCHIVE}Crick_storage/
LOVELACE=${GENOMEDARCHIVE}Lovelace_decoding/
MENDEL=${GENOMEDARCHIVE}Mendel_annotating/

# Set N of samples
nsamples=4
##___________________________##
# exit script if there is no arguments
#: ${1?"$usage"}

while getopts ':h:n:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    n) nsamples=$OPTARG
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

##################################################################
# 1) Get list of samples to run bcbio ordered by date (older to newer)
##################################################################
IFS=$' ' read -r -a sampleorder <<< $(find $(ls -dtr ${LOVELACE}Analysis/*) -maxdepth 0 -type d  \! -exec test -d '{}/bcbio' \; -print | grep -v "${LOVELACE}Analysis/$")

# ask user if want to change the order of running
echo
echo "--------------------------------------------------------"
echo "                   RUNNING ORDER                        "
echo "--------------------------------------------------------"

select sample in ${sampleorder[@]##*/} exit; do
    echo "Sample ${sample##*/} moved to the beginning of the list"

    # check for older folder, and 
    touch -t $(echo $(find ${LOVELACE}Analysis/ -maxdepth 1 -type d -printf '%Ty%Tm%Td%TH%TM.01\n' | sort | head -n1)-10000 | bc) ${LOVELACE}Analysis/${sample}

    case $sample in
      exit) 
        echo "Proceding to bcbio-nextgen analysis"
        break
        ;;		
    esac
done

# update sampleorder with user input
IFS=$' ' read -r -a sampleorder <<< $(find $(ls -dtr ${LOVELACE}Analysis/*) -maxdepth 0 -type d  \! -exec test -d '{}/bcbio' \; -print | grep -v "${LOVELACE}Analysis/$")

##################################################################
# 2) Run bcbio for multiple samples
##################################################################
# list all samples in ${LOVELACE}Analysis/
allsamples=($(ls ${LOVELACE}Analysis/))
# get samples to run, using user order
samplenames=(${sampleorder[@]##*/})

while [ ${#allsamples[@]} -ne 0  ]; do
    printf '%s\n' ${samplenames[@]} | xargs -n1 -P${nsamples} -I {} sh -c "
    echo 'Starting bcbio-nextgen analysis for {}'
    ${LOVELACE}bin/run_bcbio.sh -s {}
    echo 'Finishing bcbio-nextgen for {}'"

    for finishedsample in ${LOVELACE}Finished/*; do
      delete=(${finishedsample##*/})
      samplenames=(${samplenames[@]/$delete})
    done

    allsamples=($(ls ${LOVELACE}Analysis/))
    addedsamples=($(echo ${allsamples[@]} ${samplenames[@]} | tr ' ' '\n' | sort | uniq -u))

    for newsample in ${addedsamples[@]}; do
      echo ${newsample}
      samplenames+=(${newsample})
    done

done

