#!/usr/bin/env bash
# exit when any command fails
set -e
set -o pipefail

#######################################################################################
#										      #	
#                                                                                     #
#                            Downloading and extracting                               #
#										      #
#                                                                                     #
#######################################################################################

###################################
# Description
###################################

#download_extract.sh - download the batch using the formulation file (p.e. fromana) and extract the samples of the batch. The tmp files for the pairs will be also created in the Lovelace_decoding/Raw

###################################
#   HELP function
###################################
usage="$(basename "$0") [-h] [-f <formulation file with information: link, plataform and samples>] [-d <download tar file>] [-p <path>]
       -- program to download and extract  --
where:
    -h    show this help text
    -f    formulation file
    -d    [default: YES] YES to download the file from ana, or stats.download file name to star the analyze from extraction step
    -p    [default:/genomedarchive/] path before the Lovelace and Crick folders"

##__________ SETUP __________##
if [ -d "/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/" ]; then
  path="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
else
  path="/genomedarchive/"
fi

: ${1?"$usage"}

while getopts ':h:f:d:p:' option; do
    case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) inputfile=$OPTARG
       ;;
    d) downloadfile=$OPTARG
       ;;
    p) path=$OPTARG
       ;;
    :) printf "missing argument for -%\n" "$OPTARG" >&2
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

##__________ SET PRIMARY DIRECTORIES  __________##
path_crick=${path}Crick_storage/
path_love=${path}Lovelace_decoding/
path_rosa=${path}Rosalind_resolution/

##__________ CHECK IF ALL IS READY TO START THE PREPROCESSING ANALYSIS __________##
if [ -e ${path_crick}Source_control/gene_panel/${inputfile} ]; then 
  printf -- "Starting reading the formulation file: ${inputfile} ...\n"
else
  printf -- '\033[31m ERROR: Input file is missing \033[0m\n';
  exit 1
fi

##__________ STARTING THE DOWNLOAD AND EXTRACTION__________##
batch_download=$(head -n 1 ${path_crick}Source_control/gene_panel/${inputfile})
plataform=$(cat ${path_crick}Source_control/gene_panel/${inputfile} | head -n2 | tail -n 1)
file=$(echo ${batch_download##*/}) ##fom var, greedy front trim ##, matches anything *, until the last /
cd ${path_crick}Raw/Downloads/
if [[ $downloadfile =~ ^YES ]]; then
  printf -- "Starting the download of: ${file}...\n";
  wget -c ${batch_download} -S -r --show-progress -o stats.download
  if grep "FINISHED" stats.download; then 
    extract_name="$(grep 'Saving to' stats.download | cut -d ':' -f2 | sed -e "s/ ‘//" -e "s/’//" | uniq)"
    if [[ ${extract_name} == */* ]]; then
      mv ${extract_name} ${path_crick}Raw/Downloads/.
      folder_rm=$(echo ${extract_name} | cut -d "/" -f1)
      rm -r ${folder_rm}
      batch_name=$(echo ${extract_name##*/} | cut -d "." -f1)
      extract_name=$(echo ${extract_name##*/})
    else
      batch_name="$(grep 'Saving to' stats.download | cut -d ':' -f2 | cut -d "." -f1 | sed  "s/ ‘//g" | uniq)"
    fi
    mv stats.download ${batch_name}.stats.download
    printf -- '\033[32m SUCCESS: Download finished! \033[0m\n';
  else ##improve this, give the last files created: ls -taF | grep "tar"
    aux_var="$(head -n1 stats.download | sed "s/.*https//")"
    aux_download=$(echo ${aux_var##*/})
    printf -- "\033[33m WARNING: This batch was already download with the name of:${aux_download}  \033[0m\n";
  fi
else 
  printf -- "Passing the dowload step, for: ${downloadfile} \n" 
  extract_name="$(grep 'Saving to' ${downloadfile}  | cut -d ':' -f2 | sed -e "s/ ‘//" -e "s/’//" | uniq)"
  if [[ ${extract_name} == */* ]]; then
    mv ${extract_name} ${path_crick}Raw/Downloads/.
    folder_rm=$(echo ${extract_name} | cut -d "/" -f1)
    rm -r ${folder_rm}
    batch_name=$(echo ${extract_name##*/} | cut -d "." -f1)
    extract_name=$(echo ${extract_name##*/})
  else
    batch_name="$(grep 'Saving to' ${downloadfile}  | cut -d ':' -f2 | cut -d "." -f1 | sed  "s/ ‘//g" | uniq)"
  fi
fi
printf -- "Starting the extraction of: ${extract_name}...\n";
tar xvf ${extract_name} --skip-old-files -C ${path_love}Raw/. &> ${path_love}Raw/${batch_name}.stats.extract
if grep "skipping existing file" ${path_love}Raw/${batch_name}.stats.extract; then 
  printf -- "\033[32m SUCCESS: Extract finished but ... one warning detected: \033[0m\n";
  printf -- "\033[33m WARNING: Some files had already be extracted. The analysis will continue (possible some files will be replace)!  \033[0m\n";
else
  printf -- "\033[32m SUCCESS: Extract finished to ${path_love}Raw/${batch_name} \033[0m\n";
fi
echo "Creating samples folders, generating log reports and pairing the fasta files by samples ..."
cd ${path_love}Raw/.
folder_name="$(head -n1 ${batch_name}.stats.extract | cut -d "/" -f1)"
grep -v "skipping existing file" ${batch_name}.stats.extract | grep "fq.gz" | sort > ${folder_name}.pairs.tmp
grep -v "skipping existing file" ${batch_name}.stats.extract | grep "fq.gz" | sort > ${folder_name}.pairs.list_tmp
path_report=$(grep "eport" ${batch_name}.stats.extract | head -n1)
mkdir -p ${batch_name}_reports
cp ${path_love}Raw/${path_report}* ${batch_name}_reports/.

echo "${folder_name}.pairs.tmp $batch_name $plataform" > ${path_love}Raw/download.tmp
