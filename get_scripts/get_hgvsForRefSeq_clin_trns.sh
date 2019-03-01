#!/usr/bin/env bash
set -e 
set -o pipefail

echo "STEP1: Getting clinical RefSeq mRNA transcripts from RefSeqGRCh37_clinical_hdr_sort.bed.gz"
less RefSeqGRCh37_clinical_hdr_sort.bed.gz | cut -f10 | grep -v "refSeq_mRNA" | sort | uniq >> nm_transcripts.txt
echo "STEP1 Done."

echo "STEP2: Decompress VEP output files, if they are compressed. Retrieve RefSeq mRNA transcripts from VEP output for each chromosome..."
count=`ls chr*.vep*.cache*.homo37.vcf.gz 2>/dev/null | wc -l`
if [ $count != 0 ]
then 
gunzip chr*.vep*.cache*.homo37.vcf.gz
fi

for SAMPLE in chr1.vep95.cache95.homo37.vcf chr2.vep95.cache95.homo37.vcf chr3.vep95.cache95.homo37.vcf chr4.vep95.cache95.homo37.vcf chr5.vep95.cache95.homo37.vcf chr6.vep95.cache95.homo37.vcf chr7.vep95
.cache95.homo37.vcf chr8.vep95.cache95.homo37.vcf chr9.vep95.cache95.homo37.vcf chr10.vep95.cache95.homo37.vcf chr11.vep95.cache95.homo37.vcf chr12.vep95.cache95.homo37.vcf chr13.vep95.cache95.homo37.vcf 
chr14.vep95.cache95.homo37.vcf chr15.vep95.cache95.homo37.vcf chr16.vep95.cache95.homo37.vcf chr17.vep95.cache95.homo37.vcf chr18.vep95.cache95.homo37.vcf chr19.vep95.cache95.homo37.vcf chr20.vep95.cache9
5.homo37.vcf chr21.vep95.cache95.homo37.vcf chr22.vep95.cache95.homo37.vcf chrMT.vep95.cache95.homo37.vcf chrY.vep95.cache95.homo37.vcf chrX.vep95.cache95.homo37.vcf

do

time fgrep -f nm_transcripts.txt ${SAMPLE} >> vepout_nm_clin.tmp.txt

done
echo "STEP2 Done."

##if necessary:
echo "STEP3: Compress VEP outputs..."
for file in chr*vcf; do bgzip $file; done
echo "STEP3 Done."

#Split the output:
echo "STEP4: Split VEP output with clinical transcripts..."
split --line-bytes=10GB vepout_nm_clin.tmp.txt vepout_nm_clin.tmp.split
echo "STEP4 Done."

#Get the HGVS nomenclature for transcripts:
echo "STEP5: Get HGVS nomenclature for VEP output with clinical transcripts..."
find vepout*split* | xargs -n1 -P14 -I {} python get_hgvsnomenclature.py -vcf {} -outname {}
echo "STEP5 Done."

#Merge and sort files
echo "STEP6: Merge and Sort files..."
for file in vepout*.vcf; do more $file >> hgvsnomenclature_merged.tmp.txt; done 
for file in vepout*.warning; do more $file >> grch37.clin.hgvs_dbsnp_v0.1.warning.vcf; done

sort -k1,1 -k2,2n -T . hgvsnomenclature_merged.tmp.txt > grch37.clin.hgvs_dbsnp_v0.1.vcf
echo "STEP6 Done."

#Add header
echo "STEP7: Add header, BGZipped and Indexed output file..."
sed -i "1i\#Chro\tPos\tRef\tAlt\tRefSeqmRNA\tHGVSc\tRefSeqProtein\tHGVSp\trs_id" grch37.clin.hgvs_dbsnp_v0.1.vcf

#Zip and index files
bgzip grch37.clin.hgvs_dbsnp_v0.1.vcf
bgzip grch37.clin.hgvs_dbsnp_v0.1.warning.vcf
tabix -p vcf grch37.clin.hgvs_dbsnp_v0.1.vcf.gz
echo "STEP7 Done."

#Remove files
rm *tmp*
echo "Finished with SUCCESS!! Bye."
