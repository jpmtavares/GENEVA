#!/usr/bin/env bash

#######################################################################################
#                                                                                     #
#                 Merge vcf files and split them by chromosome                        #
#                                                                                     #
#######################################################################################

############################################################
#   SETUP
############################################################

# Provide a variable with the location of vcf files
if [ $# -eq 0 ] #check if there are any argument provided
  then
    vcfPath=$(pwd)
  else
    vcfPath="$1"
fi

############################################################
#   1) merge vcf files with vcf-merge tool
############################################################
echo "STEP 1:"
echo "Merging the vcf files ..."

# Get the running date of script
todaydate="$(date +%Y%m%d)"

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
    bcftools filter ${vcfPath}/todas_${todaydate}.vcf.gz -r ${CHROM} --threads 15 | bgzip -c > ${vcfPath}/chr${CHROM}_${todaydate}.vcf.gz;
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
find ${vcfPath}/chr*_*vcf.gz | xargs -n1 -P24 -I {} python get_inHouseFreq.py -vcf {} -outname {}tmp.freq

#Create a unique file with the frequences
cat ${vcfPath}/*tmp.freq >> ${vcfPath}/in_HouseInput_${todaydate}.txt

#Remove all temporary files
rm chr*gz*
rm *tmp.freq

echo "Finish!" 
