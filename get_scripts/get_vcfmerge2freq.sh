#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#           Merge samples VCF files and get inHouse frequency table                   #
#                                                                                     #
#######################################################################################

############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h] [-v <path to VCF files>] [-d <path to output>] [-o <output name>] [-p <number of process to run>]

       -- program to merge VCF samples, and calculate allele frequencies among those samples --

where:
    -h    show this help text
    -v    [default: current directory] set directory where VCF samples are stored
    -d    [default: same as <vcf>] set directory where allele frequency table is written
    -o    [default: inHouse_freq_] set output name for allele frequency file
    -p    [default: 24] set number of process to run at the same time"

##__________ SETUP __________##
# Get the running date of script
todaydate="$(date +%Y%m%d)"
# Set current directory as default
vcfPath=$(pwd)
# Set output directory the same as <vcf>
freqPath=${vcfPath}
# Set outputname
outputname="inHouse_freq_"
# Set number of process
processes=24

while getopts ':h:v:d:o:p:' option; do
  case "$option" in
    h) echo "$usage"
       exit
       ;;
    v) vcfPath=$OPTARG
       ;;
    d) freqPath=$OPTARG
       ;;
    o) outputname=$OPTARG
       ;;
    p) processes=$OPTARG
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
#   1) merge vcf files with vcf-merge tool
############################################################
echo "STEP 1:"
echo "Merging the vcf files ..."

echo "number of samples: "
#Count number of sample
ls ${vcfPath}/[0-9]*.vcf.gz | wc -l

# Run vcf-merge and compress file with bgzip
vcf-merge $(readlink -f ${vcfPath}/[0-9]*.vcf.gz) | bgzip -c > ${vcfPath}/todas_${todaydate}.vcf.gz
# Create tabix index file
tabix -p vcf ${vcfPath}/todas_${todaydate}.vcf.gz

echo "All the vcf files were merged"

############################################################
#   2) split merged vcf file by chromosome
############################################################
echo "STEP 2:"
echo "Split the files by chromossomes ..."

#Only the set of chromossomes were used, e.g. the MT was excluded
for CHROM in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y;
do
    bcftools filter ${vcfPath}/todas_${todaydate}.vcf.gz -r ${CHROM} | bgzip -c > ${vcfPath}/chr${CHROM}_${todaydate}.vcf.gz;
    # Create tabix index file
    tabix -p vcf ${vcfPath}/chr${CHROM}_${todaydate}.vcf.gz;
done

echo "Done this step"

############################################################
#   3) get the frequences (in house info)
############################################################
echo "STEP 3:"
echo "Construct the table with the number of samples, number of homozigotics and the allelic frequences ..."
#Get the number os samples, the homozigotics and the allelic frequences
#The number of processes runing at the same time is set to 24 (-P option)
find ${vcfPath}/chr*_*vcf.gz | xargs -n1 -P${processes} -I {} get_inHouseFreq.py -vcf {} -outname {}tmp.freq

#Create a unique file with the frequences
cat ${vcfPath}/*tmp.freq >> ${vcfPath}/in_HouseInput_${todaydate}.txt

#Sort file by chromosome position
sort -k1,1 -k2,2n ${vcfPath}/in_HouseInput_${todaydate}.txt > ${freqPath}/${outputname}${todaydate}.txt

#Remove all temporary files
rm ${vcfPath}/chr*_*vcf.gz*
rm ${vcfPath}/*tmp.freq
rm ${vcfPath}/in_HouseInput_${todaydate}.txt

echo "Finish!" 
