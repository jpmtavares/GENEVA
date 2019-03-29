#!/usr/bin/env bash
#set -e 
#set -o pipefail

#######################################################################################
#                                                                                     #
#                      Get HGVS nomenclature for RefSeq file                          #
#                                                                                     #
#######################################################################################
############################################################
#   HELP function
############################################################

usage="$(basename "$0") [-h] [-r <RefSeq file.bed.gz>] [-c <path with output chromosomes from VEP>]
       -- script that gets HGVS nomenclature for a list of RefSeq_mRNA transcripts --

where:
    -h    show this help text
    -f    [default: {LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*.txt] RefSeq file
    -c    [default: {CRICK}Annotation/Variants/VEP/] path with output chromosomes from VEP
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
##___________________________##
# Set RefSeq file
RefSeq=$(ls ${LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*.txt)

# Set VEP output chromosomes path
chromosomes=${CRICK}Annotation/Variants/VEP/

##___________________________##
# exit script if there is no arguments
#: ${1?"$usage"}

while getopts ':h' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    f) RefSeq=$OPTARG
       ;;
    c) chromosomes=$OPTARG
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


echo "_______________________________________________________________"
echo "STEP1: Getting clinical RefSeq mRNA transcripts from $(basename ${RefSeq})"
echo "_______________________________________________________________"
echo
#look for column called "refSeq_mRNA_noVersion" and print it to file 
awk -F'\t' -v c="refSeq_mRNA_noVersion" 'NR==1{for (i=1; i<=NF; i++) if ($i==c){p=i; break}; next} {print $p}' ${RefSeq} | sort | uniq >> ${MENDEL}nm_transcripts.txt

echo "Reading $(basename ${RefSeq}) and get unique transcripts in 'refSeq_mRNA_noVersion' column"

echo "_______________________________________________________________"
echo "STEP1 Done."
echo "_______________________________________________________________"

echo "_______________________________________________________________"
echo "STEP2: Uncompress VEP output files, if they are compressed."
echo "Retrieve RefSeq mRNA transcripts from VEP output for each chromosome..."
echo "_______________________________________________________________"
echo

for CHRformat in ${chromosomes}chr*.vep*.cache*.homo37.*; do
  chrname=$(basename ${CHRformat})
  echo "[" "$(date '+%Y-%m-%d %H:%M:%S' )" "]" "Get HGVS nomenclature for ${chrname%%.*}"

  # Uncompressed chromosome files if necessary
  if [ "${CHRformat: -7}" == ".vcf.gz" ]; then
    echo "                       " "Uncompress $(basename ${CHRformat})"
    gunzip ${CHRformat}
    chr=${CHRformat::-3}
    echo "                       " "$(basename ${CHRformat}) was uncompressed and save in $(basename ${chr})"
  else
    echo "                       " "$(basename ${CHRformat}) is already uncompressed."
    chr=${CHRformat}
  fi

  # Grep clinical transcripts in chromosome files
  echo "                       " "Look for clinical RefSeq mRNA transcripts in $(basename ${chr})"
  fgrep -f ${MENDEL}nm_transcripts.txt ${chr} > ${MENDEL}vepout_${chrname%%.*}.tmp.txt

  # Split VEP output
  split --line-bytes=500MB ${MENDEL}vepout_${chrname%%.*}.tmp.txt ${MENDEL}vepout_${chrname%%.*}.tmp.split

  # Run get_HGVSnomenclature.py
  echo "                       " "Run get_HGVSnomenclature.py"
  find ${MENDEL}vepout*split* | xargs -n1 -P14 -I {} python ${MENDEL}bin/get_HGVSnomenclature.py -vcf {} -outname {}

  # Merge and Sort files
  echo "                       " "Merge and Sort files"
  for file in ${MENDEL}vepout*.vcf; do
    more $file >> ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.tmp.txt;
  done

  for file in ${MENDEL}vepout*.warnings; do
    more $file >> ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.warnings.vcf; 
  done

  sort -k1,1 -k2,2n -T . ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.tmp.txt > ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.vcf

  # Zip and Index files
  echo "                       " "Zip and Index files"
  bgzip ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.vcf
  bgzip ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.warnings.vcf
  tabix -p vcf ${MENDEL}grch37.clin.hgvs_dbsnp.${chrname%%.*}.vcf.gz

  #Remove files
  rm ${MENDEL}*tmp*
  
  echo "                       " "$(basename ${chr}) is Done."
  echo

done
echo "_______________________________________________________________"
echo "STEP2 Done."
echo "_______________________________________________________________"
echo

echo "_______________________________________________________________"
echo "STEP3: Compress VEP outputs..."
echo "_______________________________________________________________"
echo

for file in ${chromosomes}chr*vcf; do 
echo "                       " "Compress $(basename $file)"
bgzip $file
echo "                       " "$(basename $file) compressed and save in $(basename $file).gz" 
echo
done

echo "_______________________________________________________________"
echo "STEP3 Done."
echo "_______________________________________________________________"

#Remove files
rm ${MENDEL}nm_transcripts.txt

echo "Finished with SUCCESS!! Bye."
