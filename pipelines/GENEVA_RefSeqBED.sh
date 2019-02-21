#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#     RefSeq       #
#                                                                                     #
#######################################################################################

# PRINT HEADER AND RUN COMMAND ON THE BODY {{{
# use it in pipelines, e.g. ps | body grep somepattern
body() {
    IFS= read -r header
    printf '%s\n' "$header"
    "$@"
}
#}}}


############################################################
#   HELP function
############################################################
usage="$(basename "$0") [-h]
       -- pipeline to create RefSeq annotation for clinical transcripts --
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
RefSeq=$(basename https://raw.githubusercontent.com/jpmtavares/GENEVA/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt)
create_RefSeqBED=$(basename https://raw.githubusercontent.com/jpmtavares/GENEVA/master/create_scripts/create_RefSeqBED.R)
genesCoordinates=$(basename https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_genesCoordinates.py)

############################################################
#   1) Download BED files from UCSC Table Browser
############################################################
# get exons.bed from UCSC Table Browser
wget --progress=dot \
'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&'\
'hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_compressType=none&outBed=1&'\
'fbQual=exon&fbExonBases=0&submit=submit&hgta_doGetBed=1' -O "exons.bed"

# get introns.bed from UCSC Table Browser
wget --progress=dot \
'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&'\
'hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_compressType=none&outBed=1&'\
'fbQual=intron&fbIntronBases=0&submit=submit&hgta_doGetBed=1' -O "introns.bed"

############################################################
#   2) Download create_script and annotation file
############################################################
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/create_scripts/create_RefSeqBED.R
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_genesCoordinates.py

# Make it executables
chmod +x ${create_RefSeqBED}
chmod +x ${genesCoordinates}

############################################################
#   3) Run create_RefSeq.R
############################################################
./create_RefSeqBED.R --exons=="exons.bed" --introns=="introns.bed" --RefSeq==${RefSeq}
#______________________________
# complete BED file - uncomment if needed
#______________________________
#mkdir RefSeq_annotation/complete/
# 3.1 sort
#less RefSeq_annotation/RefSeqGRCh37_hdr.bed | body sort -k1,1 -k2,2n - > RefSeq_annotation/complete/RefSeqGRCh37_hdr_sort.bed

# 3.2 coverage BED
#awk 'NF{NF-=1};1' RefSeq_annotation/complete/RefSeqGRCh37_hdr_sort.bed | uniq - | awk '{ print $1"\t"$2"\t"$3"\t"$6","$4","$5","$8","$9","$10}' - | sed -r 's/\s+/\t/g' - > RefSeq_annotation/complete/RefSeqGRCh37_coverage.bed

# compress file and tabix
#bgzip RefSeq_annotation/complete/RefSeqGRCh37_hdr_sort.bed
#tabix -b 2 -e 3 RefSeq_annotation/complete/RefSeqGRCh37_hdr_sort.bed.gz

#______________________________
# clinical BED file
#______________________________
mkdir RefSeq_annotation/clinical/
# 3.1 sort
less RefSeq_annotation/RefSeqGRCh37_clinical_hdr.bed | body sort -k1,1 -k2,2n - > RefSeq_annotation/clinical/RefSeqGRCh37_clinical_hdr_sort.bed

# 3.2 coverage BED
awk 'NF{NF-=1};1' RefSeq_annotation/clinical/RefSeqGRCh37_clinical_hdr_sort.bed | uniq - | awk '{ print $1"\t"$2"\t"$3"\t"$6","$4","$5","$8","$9","$10}' - | sed -r 's/\s+/\t/g' - > RefSeq_annotation/clinical/RefSeqGRCh37_clinical_coverage.bed

# compress file and tabix
bgzip RefSeq_annotation/clinical/RefSeqGRCh37_clinical_hdr_sort.bed
tabix -b 2 -e 3 RefSeq_annotation/clinical/RefSeqGRCh37_clinical_hdr_sort.bed.gz

############################################################
#   4) Run get_genesCoordinates.py for complete (uncomment if needed) and clinical BED files
############################################################
#./get_genesCoordinates.py --refseq RefSeq_annotation/complete/RefSeqGRCh37_hdr_sort.bed.gz --out_file RefSeq_annotation/complete/RefSeqGRCh37_coordinates.txt

./get_genesCoordinates.py --refseq RefSeq_annotation/clinical/RefSeqGRCh37_clinical_hdr_sort.bed.gz --out_file RefSeq_annotation/clinical/RefSeqGRCh37_clinical_coordinates.txt

############################################################
#   5) Remove intermediate files
############################################################
rm exons.bed
rm introns.bed
rm RefSeq_annotation/RefSeqGRCh37.bed
rm RefSeq_annotation/RefSeqGRCh37_clinical.bed
rm ${create_RefSeqBED}
rm ${genesCoordinates}