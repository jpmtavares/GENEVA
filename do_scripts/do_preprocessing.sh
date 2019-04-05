#!/usr/bin/env bash
# exit when any command fails
set -e
set -o pipefail

#######################################################################################
#										      #	
#                                                                                     #
#                            Organizing and preprocessing                             #
#										      #
#                                                                                     #
#######################################################################################

###################################
# Description
###################################

#samples_preprocessing.sh - given a sorted file (p.e. pairs.tmp files) with the samples to organize (create a fold with the sample name in Lovelace_decoding/Analysis) and preprocessing with fastp 
    #Used scripts: get_sample_name.py and fastp

###################################
#   HELP function
###################################
usage="$(basename "$0") [-h] [-f <input file with the samples to analyze>] [-b <batch name>] [-p <plataform>] [-t <number of process to run>] 
       -- program to preprocessing the reads (after organize the samples on input file)  --
where:
    -h    show this help text
    -f    tmp file with the pairs to analyze
    -b    batch name
    -p    plataform used
    -t    [default: 5] set number of process to run - tool: fastp"

##__________ SETUP __________##
if [ -d "/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/" ]; then
  path="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
else
  path="/genomedarchive/"
fi
#path="/mnt/data/Genomed_server/"
## Set number of threads, process
processes=5

: ${1?"$usage"}

while getopts ':h:f:b:p:t:' option; do
    case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) inputfile=$OPTARG
       ;;
    b) batchname=$OPTARG
       ;;
    p) plataform=$OPTARG
       ;;
    t) processes=$OPTARG
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
if [ -e ${path_love}Raw/${inputfile} ]; then 
  printf -- "Starting reading the pairs and creating the folders...\n"
  folder_name=$( echo ${path_love}Raw/${inputfile} | cut -d "." -f1)
  total_pairs=$( more ${path_love}Raw/${inputfile} | wc -l)
  printf -- "A total number of pairs will be analyzed: "$((total_pairs/2))"\n" 
else
  printf -- '\033[31m ERROR: Input file is missing \033[0m\n';
  exit 1
fi
cd ${path_love}Raw/
while read -r diret; do
  aux=$(dirname $diret)
  count=$(ls $aux/*gz | wc -l)
  if [[ ${count} -ge 4 ]]; then
    printf -- "\033[33m WARNING: More than one pair for the same sample! \033[0m\n"; 
    cd ${aux}
    ls *1.fq.gz > tmp.files
    for file in *1.fq.gz; do zcat $file | echo $((`wc -l`/4)) >> tmp.listReads; done
    moreReads=$(python ${path_love}bin/getPairsWithMoreReads.py -list tmp.listReads)
    getPair=$(sed -n "${moreReads}p" tmp.files | sed s'/1.fq.gz//')
    rmPairs=$(sed -n "${moreReads}p" tmp.files | cut -d "_" -f1,2)
    grep "${getPair}" ${path_love}Raw/${inputfile} > tmp.files
    sed -i "/${rmPairs}/d" ${path_love}Raw/${inputfile} 
    listtmp_batch=$(echo ${inputfile/tmp/list_tmp})
    sed -i "/${rmPairs}/d" ${path_love}Raw/${listtmp_batch}
    cat tmp.files >> ${path_love}Raw/${inputfile}
    cat tmp.files >> ${path_love}Raw/${listtmp_batch}
    rm tmp*
  fi
  cd ${path_love}Raw/
done < ${inputfile}

##__________ CHECK AND GET THE SAMPLE NAME __________##

cd ${path_love}Raw/

while read -r first_pair; do
  read -r second_pair
  var=$(python ${path_love}bin/get_sample_name.py -f ${first_pair} -r ${second_pair} -p ${plataform} -b ${batchname} -a ${path_crick})
  if [[ $var =~ ^ERROR ]]; then
    sed -i '$ d' ${path_crick}Source_control/sample_batch.info
    printf -- "\033[31m ERROR: ${batchname} need to be checked! \033[0m\n";
    printf -- "\033[31m ${var} \033[0m\n";
    printf -- "\033[33m WARNING: This sample will not be preprocessed, starting preprocessing for the next sample! \033[0m\n";
    continue
  fi
  IFS="," read -ra elements <<< "$var"
  check_sample="${elements[6]}"
  sample_name="${elements[0]}"
  if [[ $check_sample =~ ^samplexist ]]; then
    printf -- "\033[31m ERROR: ${sample_name} need to be checked! \033[0m\n";
    printf -- "\033[31m STOPPING THE ANALYSIS FOR: ${sample_name} This sample already exists in sample_batch.info. \033[0m\n";
    rm log_
    continue
  fi
  if [[ $check_sample =~ ^checkbigger ]]; then
    printf -- "\033[33m WARNING: This ${sample_name} is duplicate, check for the pairs with more coverage. \033[0m\n";
    printf -- "\033[31m STOPPING THE ANALYSIS FOR: ${sample_name} This sample already exists in sample_batch.info, in the same batch. Please check. \033[0m\n";
    rm log_
    continue
  fi
  out_f="${elements[1]}"
  out_r="${elements[2]}"
  pair1="${elements[3]}"
  pair2="${elements[4]}"
  raw_name="${elements[5]}"
  if [[ $sample_name == "Empty" ]]; then
     printf -- "The sample name has a format different from usual, please insert the sample name in the following format (NNNNN-LL):\n"
     read -p  "New sample name: " </dev/tty
     sample_name="$REPLY"
     out_f=${sample_name}"_1.fastp.fq.gz"
     out_r=${sample_name}"_2.fastp.fq.gz"
     sed -i "s/Empty/${sample_name}/g" ${path_crick}Source_control/sample_batch.info
     printf -- "\033[33m WARNING: Please add this new format: ${sample_name} \033[0m\n";
     touch ${path_love}to_do/add_new_format_sample_name.txt
     echo ${first_pair} >> ${path_love}to_do/new_format_sampleName.add
  fi
  echo "$sample_name" >> ${path_crick}Source_control/${batchname}.batch_samplesList
  log_name="log_"$sample_name
  if [ -e ${path_love}log/${log_name} ]; then 
     printf -- "\033[33m WARNING: The log for ${sample_name} already exist! The analysis will continue, some files will be replaced by new files from this round! \033[0m\n";
     echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Warning: Preparing the reads for the preprocessing" >> ${path_love}log/$log_name
     echo "                        ==> First pair: $out_f" >> ${path_love}log/$log_name
     echo "                        ==> First pair: $out_r" >> ${path_love}log/$log_name
  else
  mv log_ $log_name
  echo "                        ==> First pair: $out_f" >> $log_name
  echo "                        ==> First pair: $out_r" >> $log_name
  mv $log_name ${path_love}log/.
  fi
  var=$(echo ${first_pair} | sed -e "s/\/${pair1}//g")
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Starting the preprocessing for ${sample_name}"
  echo "                        ==> Pairs are correct, so the files will be moved."
  directory=$(dirname ${first_pair})
  mkdir -p ${sample_name}
  cp ${path_love}Raw/${directory}/* ${sample_name}/.  
  rm -r ${path_love}Raw/${directory}
  echo "                        ==> All organized, starting the fastp tool."
  mkdir -p ${path_love}Analysis/${sample_name}
  mkdir -p ${path_love}Analysis/${sample_name}/fastp
  cd ${path_love}Analysis/${sample_name}/fastp
  log_fastp="fastp_log_"${sample_name}
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Starting the preprocessing with FASTP tool:" >> ${path_love}log/${log_name}
  fastp -i ${path_love}Raw/$sample_name/${pair1} -o ${out_f} -I ${path_love}Raw/$sample_name/${pair2} -O ${out_r} -w ${processes} -z 2 &>> ${log_fastp}
  echo "                        ==> Fastp tool finished, checking if an error occur!"
  if tail ${log_fastp} | grep "time used" ; then
    printf -- '\033[32m                        SUCCESS: Fastp finished! \033[0m\n';
    read1_b=$(grep -A1 "Read1 before filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read1_a=$(grep -A1 "Read1 after filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read2_b=$(grep -A1 "Read2 before filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read2_a=$(grep -A1 "Read2 aftering filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    version=$(grep "fastp v" $log_fastp | cut -d "," -f1)
    pair1_perc=$(echo "scale=4; $read1_a/$read1_b*100" | bc)
    pair2_perc=$(echo "scale=4; $read2_a/$read2_b*100" | bc)
    echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Results from fastp:" >> ${path_love}log/${log_name}
    echo "[   Software Version  ]" ${version} >> ${path_love}log/${log_name}
    echo "                        Total number of reads before preprocessing:" >> ${path_love}log/${log_name}
    echo "                        ==> First pair: $read1_b reads" >> ${path_love}log/${log_name}
    echo "                        ==> Second pair: $read2_b reads" >> ${path_love}log/${log_name}
    echo "                        Total numer of reads after preprocessing:" >> ${path_love}log/${log_name}
    echo "                        ==> First pair: $read1_a reads($pair1_perc %)" >> ${path_love}log/${log_name}
    echo "                        ==> Second pair: $read2_a reads ($pair2_perc %)" >> ${path_love}log/${log_name}
  else
    printf -- "\033[31m                        ERROR: Check log from fastp for this sample: ${sample_name} \033[0m\n";
    echo "[        ERROR!!      ]" "Check log from fastp" >> ${path_love}log/${log_name}
    mv ${log_fastp} ${path_love}log/.
    rm -r ${path_love}Analysis/${sample_name}/
    continue
  fi
  cd ${path_love}Raw/.
  mv -n ${sample_name} ${path_crick}Raw/2019/.
  if [ -d "${path_love}Analysis/${sample_name}" ]; then
    sed -i "/${pair1}/d" ${folder_name}.pairs.list_tmp
    sed -i "/${pair2}/d" ${folder_name}.pairs.list_tmp
  fi
done < ${folder_name}.pairs.tmp
cd ${path_love}Raw/
lines=$(more ${folder_name}.pairs.list_tmp | wc -l)
if [ ${lines} -eq "0" ]; then # nao funciona if [ -s ${folder_name}.pairs.list_tmp ]; then
  printf -- '\033[32m                        SUCCESS: All samples in this batch finished the preprocessing with sucess, so all the tmp or unnecessary files will be removed! \033[0m\n';
  mv *_reports/* ${path_crick}Raw/2019/Reports/.
  rm -r ${folder_name}
else
  printf -- "\033[33m WARNING: Some files from ${batchname} weren\'t removed, check the reason! \033[0m\n";
fi
