#!/usr/bin/env bash
# exit when any command fails
set -e
set -o pipefail

#######################################################################################
#										      #	
#                                                                                     #
#     Start the pipeline: download the batch, get samples, genes panel and reports    #
#										      #
#                                                                                     #
#######################################################################################

###################################
# Description
###################################

#GENEVA_getStarting - pipeline to download the batch, extract the samples in the batch, preprocessing the pair fastas and associating the gene panel for each samples downloaded. This scripts also organize the samples, moving them to their respective directory and modifyng modifying their name. A log file is created to report all the steps, softawres and files version, also a resume of some results  are reported in the log file. The preprocessing step is done with Fastp (default parameters) and the genes in the gene panel are verified, checked and associated to their respective sample. If all finished with sucess the samples were moved to the Crick_Storage and protected with sudo.
    #Used scripts: download_extract.sh, samples_preprocessing and get_genepanel.R  

###################################
#   HELP function
###################################
usage="$(basename "$0") [-h] [-f <formulation file with information: link, plataform and samples>] [-p <password to run the chown>] [-t <number of process to run>] 
       -- program that read a formulation file with information regarding a batch and organize the information to start the analysis  --
where:
    -h    show this help text
    -f    formulation file
    -p    password to run the chown for raw data
    -t    [default: 5] set number of process to run  - tool: fastp"

##__________ SETUP __________##
if [ -d "/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/" ]; then
  path="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
else
  path="/genomedarchive/"
fi

## Set number of threads, process
processes=5

: ${1?"$usage"}

while getopts ':h:f:t:p:' option; do
    case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) inputfile=$OPTARG
       ;;
    p) pass=$OPTARG
       ;;
    t) processes=$OPTARG
       ;;
    :) printf -- "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
   \?) printf -- "illegal option: -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1
       ;;
  esac
done
shift $((OPTIND - 1))

##__________ SET PRIMARY DIRECTORIES  __________##
path_crick=${path}Crick_storage/
path_love=${path}Lovelace_decoding/
path_rosa=${path}Rosalind_resolution/

##_______________CHECK IF ALL IS READY TO START THE PREPROCESSING ANALYSIS_______________##
if [ -e ${path_crick}Source_control/gene_panel/${inputfile} ]; then 
  printf -- "Starting the analysis for this batch ...\n"
else
  printf -- '\033[31m ERROR: Input file is missing \033[0m\n';
  exit 1
fi
_=$(command -v ${path_love}bin/get_genepanel.R);
if [ "$?" != "0" ]; then
  printf -- '\033[31m ERROR: Missing the script: get_genepanel.R, check this error! \033[0m\n';
  exit 1
fi
_=$(command -v ${path_love}bin/get_sample_name.py);
if [ "$?" != "0" ]; then
  printf -- '\033[31m ERROR: Missing the script: get_sample_name.py, check this error! \033[0m\n';
  exit 1
fi
_=$(command -v ${path_love}bin/download_extract.sh);
if [ "$?" != "0" ]; then
  printf -- '\033[31m ERROR: Missing the script: download_extract.sh, check this error! \033[0m\n';
  exit 1
fi

##__________ STARTING THE DOWNLOAD AND EXTRACTION__________##
${path_love}bin/download_extract.sh -f ${inputfile}
line=$(head -n 1 ${path_love}Raw/download.tmp)
IFS=' ' read -ra elements <<< "${line}"
pairs_file="${elements[0]}"
batch_name="${elements[1]}"
plataform="${elements[2]}"
rm ${path_love}Raw/download.tmp
##__________ STARTING THE PREPROCESSING__________##
${path_love}bin/do_preprocessing.sh -f ${pairs_file} -b ${batch_name} -p ${plataform}
mv ${path_crick}Source_control/${batch_name}.batch_samplesList ${path_crick}Source_control/${inputfile}.batch_samplesList
more  ${path_crick}Source_control/${inputfile}.batch_samplesList

##__________ GETTING THE GENES OF GENE PANEL ___________##
echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Getting the genes in gene panel for each sample in the ${inputfile}"
${path_love}bin/get_genepanel.R --fromana=="${inputfile}"

##__________ PROTECTING THE RAW DATA AND CHECK IF AN ERROR OCCURRED __________##
cd ${path_crick}Source_control/
tail -n +3 ${path_crick}Source_control/gene_panel/${inputfile} | cut -f1 > ${path_crick}Source_control/${inputfile}.genepanel_samplesList
while IFS='\n' read -r line || [[ -n "$line" ]]; do
  if ! grep -q "Starting Gene Panel" ${path_love}log/log_${line}; then 
    printf -- "\033[33m WARNING: ${line} don't have an gene panel associated! Check. \033[0m\n";
  else
    if grep -q "ERROR" ${path_love}log/log_${line}; then 
      errortype=$(grep -A1 "ERROR" ${path_love}log/log_${line} | tail -n1 | sed -zE 's/[[:space:]]*([[:space:]])/\1/g') 
      printf -- "\033[31m ERROR: ${line} \033[0m\n"
      printf -- "\033[31m       ${errortype} \033[0m\n"
    fi
  fi  
if [ -f ${path_love}to_do/${line}_genepanel.pending ]; then 
  printf -- "\033[33m WARNING: ${line} have genes pending! Go to to_do and check. \033[0m\n"; 
fi
done < *${inputfile}*batch_samplesList*

if grep -q -v -F -x -f ${inputfile}.batch_samplesList ${inputfile}.genepanel_samplesList; then
  samplesnobatch=$(grep -v -F -x -f ${inputfile}.batch_samplesList ${inputfile}.genepanel_samplesList | sed -zE 's/\n/    /g')
  printf -- "\033[31m ERROR: Sample(s) in ${inputfile} but not in batch. \033[0m\n"
  printf -- "\033[31m        ${samplesnobatch} \033[0m\n"
else
  printf -- "\033[32m SUCCESS: All samples have a gene panel associated \033[0m\n";
fi
rm ${path_crick}Source_control/${inputfile}.genepanel_samplesList
rm ${path_crick}Source_control/${inputfile}.batch_samplesList
rm ${path_crick}Source_control/gene_panel/${inputfile}

cd ${path_crick}Raw/2019/
echo "${pass}" | sudo -S chown -R root:root *-*

printf -- "\n-----------------------------------------------------------------\n"
printf -- "Finished this part of the analysis! Ready to start the bcbio part\n"
