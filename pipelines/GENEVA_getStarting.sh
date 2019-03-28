#!/usr/bin/env bash
# exit when any command fails
set -e
#set -o pipefail

#######################################################################################
#                                                                                     #
#     Start the pipeline: download the batch, get samples, genes panel and reports    #
#                                                                                     #
#######################################################################################

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h] [-f <formulation file with information: link, plataform and samples [-t <number of process to run>] [-p <password to run the chown> ]
       -- program that read a formulation file with information regarding a batch and organize the information to start the analysis  --
where:
    -h    show this help text
    -f    formulation file
    -t    [default: 5] set number of process to run at the same time - tool: fastp
    -p    password to run the chown for raw data"
##__________ SETUP __________##
## Set number of threads, process
processes=5

while getopts ':h:f:t:p:' option; do
    case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) inputfile=$OPTARG
       ;;
    t) processes=$OPTARG
       ;;
    p) pass=$OPTARG
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
if [ -d "/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/" ]; then
  path="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
else
  path="/genomedarchive/"
fi
path_crick=${path}Crick_storage/
path_love=${path}Lovelace_decoding/

echo "##__________ Reading the formulation file  __________##"
#download_samples=$(head -n 1 ${inputfile})
#file=$(echo ${download_samples##*/}) ##fom var, greedy front trim ##, matches anything *, until the last /
#plataform=$(cat ${inputfile} | head -n2 | tail -n 1)
cd ${path_crick}Raw/Downloads/

#echo "Download the samples ..."
#echo "Start download and extract files ..."
#wget -c $download_samples -S &> stats.download
batch_name="$(grep 'Saving to' stats.download | cut -d ':' -f2 | cut -d "." -f1 | sed  "s/ ‘//g" | uniq)"
extract_name="$(grep 'Saving to' stats.download | cut -d ':' -f2 | sed -e "s/ ‘//" -e "s/’//" | uniq)"
mv stats.download ${batch_name}.stats.download
echo "Extract the samples ..."
tar xvf ${extract_name} --skip-old-files -C ${path_love}Raw/. &> ${path_love}/Raw/${batch_name}.stats.extract

#batch_name="C204HW19020118_20190228_02_s5UxKR"
plataform="NovoGene"
echo "Start reading pairs and update or writing logs ..."
cd ${path_love}Raw/.
folder_name="$(head -n1 ${batch_name}.stats.extract | cut -d "/" -f1)"
grep "fq.gz" ${batch_name}.stats.extract | sort > ${folder_name}.pairs.tmp
grep "fq.gz" ${batch_name}.stats.extract | sort > ${folder_name}.pairs.list_tmp
path_report=$(grep "eport" ${batch_name}.stats.extract | head -n1)
mkdir -p ${batch_name}_reports
cp ${path_love}Raw/${path_report}* ${batch_name}_reports/.

while read -r first_pair; do
  read -r second_pair
  var=$(python ${path_love}bin/get_sample_name.py -f ${first_pair} -r ${second_pair} -p ${plataform} -b ${batch_name})
  if [[ $var =~ ^ERROR ]]; then
    echo " ----- !!!!An error!!!! -----"
    echo ${var}
    sed -i '$ d' ${path_crick}Source_control/sample_batch.info
    echo "From the batch ${batch_name} need to be checked!"
    echo " ----- !!!!An error!!!! -----"
    continue
  fi
  IFS="," read -ra elements <<< "$var"
  sample_name="${elements[0]}"
  out_f="${elements[1]}"
  out_r="${elements[2]}"
  pair1="${elements[3]}"
  pair2="${elements[4]}"
  raw_name="${elements[5]}"
  if [[ $sample_name == "Empty" ]]; then
     echo "The sample name has a format different from usual, please insert the sample name in the following format (NNNNN-LL):"
     read -p  "--> " </dev/tty
     sample_name="$REPLY"
     out_f=${sample_name}"_1.fastp.fq.gz"
     out_r=${sample_name}"_2.fastp.fq.gz"
  fi
  log_name="log_"$sample_name
  mv log_ $log_name
  echo "                        ==> First pair: $out_f" >> $log_name
  echo "                        ==> First pair: $out_r" >> $log_name
  mv $log_name ${path_love}log/.
  var=$(echo ${first_pair} | sed -e "s/\/${pair1}//g")
  echo "Pairs are correct, so the files will be moved ..."
  directory=$(dirname ${first_pair})
  mkdir -p ${sample_name}
  cp ${path_love}Raw/${directory}/* ${sample_name}/.
  rm -r ${path_love}Raw/${directory}

  echo "All organized, starting the fastp tool ... for ${sample_name}"
  mkdir -p ${path_love}Analysis/${sample_name}
  mkdir -p ${path_love}Analysis/${sample_name}/fastp
  cd ${path_love}Analysis/${sample_name}/fastp
  log_fastp="fastp_log_"${sample_name}
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Starting the preprocessing with FASTP tool:" >> ${path_love}log/${log_name}
  fastp -i ${path_love}Raw/$sample_name/${pair1} -o ${out_f} -I ${path_love}Raw/$sample_name/${pair2} -O ${out_r} -w ${processes} -z 2 &>> ${log_fastp}
  echo "Fastp tool finished, check if an error occur!"
  if tail ${log_fastp} | grep "time used"; then
    read1_b=$(grep -A1 "Read1 before filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read1_a=$(grep -A1 "Read1 after filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read2_b=$(grep -A1 "Read2 before filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    read2_a=$(grep -A1 "Read2 aftering filtering" $log_fastp | grep "total" | cut -d ":" -f2 | tr -d '[:space:]')
    version=$(grep "fastp v" $log_fastp | cut -d "," -f1)
    pair1_perc=$(echo "scale=4; $read1_a/$read1_b*100" | bc)
    pair2_perc=$(echo "scale=4; $read2_a/$read2_b*100" | bc)
    echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Result from $version:" >> ${path_love}log/${log_name}
    echo "                        Total numer of reads before preprocessing:" >> ${path_love}log/${log_name}
    echo "                        ==> First pair: $read1_b reads" >> ${path_love}log/${log_name}
    echo "                        ==> Second pair: $read2_b reads" >> ${path_love}log/${log_name}
    echo "                        Total numer of reads after preprocessing:" >> ${path_love}log/${log_name}
    echo "                        ==> First pair: $read1_a reads($pair1_perc %)" >> ${path_love}log/${log_name}
    echo "                        ==> Second pair: $read2_a reads ($pair2_perc %)" >> ${path_love}log/${log_name}
  else
    echo " ----- !!!!An error!!!! -----"
    echo "Some error occurs when running the fastp for the sample ${sample_name}, this need to be checked!"
    echo " ----- !!!!An error!!!! -----"
    mv ${log_fastp} ${path_love}log/.
    rm -r ${path_love}Analysis/${sample_name}/
    continue
  fi
  cd ${path_love}Raw/.
  echo "${sample_name}"

  mv -n ${sample_name} ${path_crick}Raw/2019/.
  if [ -d "${path_love}Analysis/${sample_name}" ]; then
    sed -i "s/.*${sample_name}//g" ${folder_name}.pairs.list_tmp
    sed -i "/^[[:space:]]*$/d" ${folder_name}.pairs.list_tmp
  fi
done < ${folder_name}.pairs.tmp

echo "Creating the gene panel for each sample ..."
Rscript ${path_love}bin/get_genepanel.R --fromana=="${inputfile}"
cd ${path_love}Raw/.
if [ -s ${folder_name}.pairs.list_tmp ]; then
  cd ${path_crick}Source_control/gene_panel 
  if find . -name '*pending' -printf 1 -quit | grep -q 1; then 
    ls *genepanel.pending >> pending_file.tmp
    while IFS='' read -r line || [[-n "$line" ]]; do
      sample_id=$(echo ${line} | cut -d "_" -f1)
      mv ${line} ${path_lovelace}Raw/.
      mv ${path_crick}Raw/2019/${sample_id} ${path_lovelace}Raw/.
    done < pending_file.tmp
    rm pending_file.tmp 
  fi
  ls *genepanel.txt >> genepanels_file.tmp
  while IFS='' read -r line || [[ -n "$line" ]]; do
    sample_id=$(echo ${line} | cut -d "_" -f1)   
    mv ${line} ${path_crick}Raw/2019/${sample_id}/.
  done < genepanels_file.tmp
  rm genepanels_file.tmp    
  mv ${path_love}Raw/${batch_name}_reports/* ${path_crick}Raw/2019/Reports/.
  rm -r ${path_love}Raw/${folder_name}*tmp
  rm -r ${path_love}Raw/${batch_name}*
else
  if find . -name '*pending' -printf 1 -quit | grep -q 1; then
    ls *genepanel.pending >> pending_file.tmp
    while IFS='' read -r line || [[ -n "$line" ]]; do
      sample_id=$(echo ${line} | cut -d "_" -f1)
      if [ -d "${path_crick}Raw/2019/${sample_id}" ]; then
        mv ${path_crick}Raw/2019/${sample_id} ${path_lovelace}Raw/.
      fi 
      mv ${line} ${path_lovelace}Raw/.
    done < pending_file.tmp
    rm pending_file.tmp
  else
    ls *genepanel.txt >> genepanels_file.tmp
    while IFS='' read -r line || [[ -n "$line" ]]; do
      sample_id=$(echo ${line} | cut -d "_" -f1)
      if [ -d "${path_crick}Raw/2019/${sample_id}" ]; then
        mv ${line} ${path_crick}Raw/2019/${sample_id}/.
      else
        mv ${line} ${path_lovelace}Raw/.
      fi
    done < genepanels_file.tmp
    rm genepanels_file.tmp
  fi
fi
cd ${path_love}Raw/
rm -r ${folder_name}
cd ${path_crick}Raw/2019/.
echo "${pass}" | sudo -S chown -R root:root *-*

echo "##__________ Extract and Preprocessing done! __________##"
echo "##__________  Good luck for the next steps!  __________##"
