#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#     Download and Run get_vcfmerge2freq.sh and get_inHouseFreq.py from GENEVA        #
#                                                                                     #
#######################################################################################

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h]

       -- pipeline to get Allele Frequencies from GENEVA --
       check https://github.com/jpmtavares/GENEVA for more information


where:
    -h    show this help text
"

while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
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

############################################################
#   SETUP
############################################################
vcfmerge2freq=$(basename https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_vcfmerge2freq.sh)
inHouseFreq=$(basename https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_inHouseFreq.py)

if [ -d "/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/" ]; then
  path="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
else
  if [ -d "/genomedarchive/" ]; then
    path="/genomedarchive/"
  else
    path="/mnt/data/Genomed_server/"
  fi
fi

############################################################
#   SET PRIMARY DIRECTORIES
############################################################
path_crick=${path}Crick_storage/
path_love=${path}Lovelace_decoding/
path_rosa=${path}Rosalind_resolution/
path_mend=${path}Mendel_annotating/
path_arch=${path}Archive/

############################################################
#   1) Create folder and symbolic link for samples (if necessary)
############################################################
cd ${path_mend}inHouse-VCF

if ! [ -d "Repository" ]; then
  mkdir Repository
  cd Repository/
  #From the Archive folder
  for file in ../../../Archive/VCF/*-gatk-haplotype*vcf.gz*; do ln -s $file .; done
  #Create info_samples to track the samples that already had a symbolic link
  for file in [0-9]*.vcf.gz; do echo $file | cut -d "-" -f1,2 | sort | uniq >> info_samples.txt; done
  #From the Archive folder, 2018 year
  for file in ../../../Archive/Analysis/2018/NovoGene/*; do echo $file >> samples2018.tmp.txt; done
  while read -r diret; do
    samplename=${diret##*/}
    if ! grep --quiet "${samplename}" info_samples.txt; then
      echo $samplename " ....... doesnt exists"
      find ${diret}/final/${samplename}/ -name '*-gatk-haplotype.vcf.gz*' >> ${samplename}.tmp
      while read -r line; do
        ln -s ${line} .
      done < ${samplename}.tmp
      rm ${samplename}.tmp
      echo ${samplename}
    fi
    cd ${path_mend}inHouse-VCF/Repository
  done < samples2018.tmp.txt
  rm samples2018.tmp.txt
  rm info_samples.txt
  for file in [0-9]*.vcf.gz; do echo $file | cut -d "-" -f1,2 | sort | uniq >> info_samples.txt; done
  #From the Archive folder, 2019 year
  #exception: if necessary overwrite an existing symlink
  ln -sf ../../../Archive/Analysis/2019/NovoGene/68277-1-/work/gatk-haplotype/68277-1--effects-annotated-filterSNP-filterINDEL.vcf.gz 68277-1--gatk-haplotype.vcf.gz
  ln -sf ../../../Archive/Analysis/2019/NovoGene/68277-1-/work/gatk-haplotype/68277-1--effects-annotated-filterSNP-filterINDEL.vcf.gz.tbi 68277-1--gatk-haplotype.vcf.gz.tbi
  rm 50927-NN-gatk-haplotype-annotated.vcf*
  rm 39069-AT-gatk-haplotype-annotated.vcf*
  for file in ../../../Archive/Analysis/2019/*/*; do echo $file >> samples2019.tmp.txt; done
  while read -r diret; do
    name_aux=${diret##*/} #because the samples name with -1-
    samplename=$(echo $name_aux | cut -d "-" -f1,2)
    folder=$(echo $diret | rev | cut -d '/' -f2 | rev)
    if [[ ${folder} != "2nd_level" ]] ; then
      if ! grep --quiet "${samplename}" info_samples.txt; then
      echo "${samplename}" >> info_samples.txt
      echo $samplename " ....... doesnt exists"
      find ${diret}/final/ -name '*-gatk-haplotype.vcf.gz*' >> ${samplename}.tmp
      while read -r line; do
        ln -s ${line} .
      done < ${samplename}.tmp
      rm ${samplename}.tmp
      echo ${samplename}
      fi
    cd ${path_mend}inHouse-VCF/Repository
    fi
  done < samples2019.tmp.txt
  rm samples2019.tmp.txt
  rm info_samples.txt
fi
cd ${path_mend}inHouse-VCF/Repository/
############################################################
#   2) Download get_scripts
############################################################
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_vcfmerge2freq.sh
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_inHouseFreq.py

# Make them executables
chmod +x ${vcfmerge2freq}
chmod +x ${inHouseFreq}

############################################################
#   3) Get the last samples from 2019
############################################################
if [ ! -f info_samples.txt ]; then 
  for file in [0-9]*.vcf.gz; do echo $file | cut -d "-" -f1,2 | sort | uniq >> info_samples.txt; done
fi
#From the Crick folder
samplesmanual=$(ls ../../../Crick_storage/Analysis/*/ | grep "manual")
for file in ../../../Crick_storage/Analysis/2019/*; do echo $file >> samplescrick.tmp.txt; done
sed -i "/manual/d" samplescrick.tmp.txt

printf -- "\033[33mWARNING: Some folders were manual obtained \033[0m\n";
printf -- "\033[33m${samplesmanual}  \033[0m\n"; 
while read -r diret; do
  name_aux=${diret##*/} #because the samples name with -1-
  samplename=$(echo $name_aux | cut -d "-" -f1,2)
  samplename=${diret##*/}
  echo $samplename
  if ! grep --quiet "${samplename}" info_samples.txt; then
    echo "${samplename}" >> info_samples.txt
    echo $samplename " ....... doesnt exists"
    find ${diret}/bcbio/*/final/ -name '*-gatk-haplotype*.vcf.gz*' >> ${samplename}.tmp
    while read -r line; do
      ln -s ${line} .
    done < ${samplename}.tmp
    rm ${samplename}.tmp
    echo ${samplename}
  cd ${path_mend}inHouse-VCF/Repository
  fi
  cd ${path_mend}inHouse-VCF/Repository
done < samplescrick.tmp.txt
rm samplescrick.tmp.txt
rm info_samples.txt
for file in [0-9]*.vcf.gz; do echo $file | cut -d "-" -f1,2 | sort | uniq >> info_samples.txt; done

############################################################
#   3) Run get_vcfmerge2freq.sh
############################################################
# Run with default parameters: present work directory will be used to get VCF files and write Allele frequency output file
./get_vcfmerge2freq.sh

############################################################
#   4) Remove folder
############################################################
cd ${path_mend}inHouse-VCF/Repository/
rm ${vcfmerge2freq}
rm ${inHouseFreq}
