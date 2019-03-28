#!/usr/bin/env bash
set -e 
set -o pipefail
#######################################################################################
#                                                                                     #
#                            Build RefSeq Annotation files                            #
#                                                                                     #
#######################################################################################
#                                                                                     #
# This script builds RefSeq annotation files needed for some downstream analysis. It  #
# starts by evaluating if there's some pending update of clinical transcripts in      # 
# ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update, since this information   #
# will be considered in several other annotation files.                               #  
#                                                                                     #
# REQUIRED:                                                                           #
#   - create_grch37.refseq_ensembl_lrg_hugo.R                                         #
#   - create_grch37.exons.refseq_ensembl_lrg_hugo.R                                   #
#   - get_genesCoordinates.py                                                         #
#                                                                                     #
# 1) create_grch37.refseq_ensembl_lrg_hugo.R                                          #
# REQUIRED:                                                                           #
#   - NCBI transcripts (downloaded from: https://www.ncbi.nlm.nih.gov/projects/       #
#genome/guide/human/index.shtml)                                                      # 
#   - Ensembl.grch37vs38_genes.txt (from: GENEVA)                                     #
#   - Ensembl BIOMART (R package)                                                     #
#   - LRG (from: ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_GRCh37.bed)              #
#   - grch37.clin.manual.refseq_ensembl.txt (from: GENEVA)                            #
#   - HUGO genes (from: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/          #
#locus_types/gene_with_protein_product.txt)                                           # 
#                                                                                     #
# 2) create_grch37.exons.refseq_ensembl_lrg_hugo.R                                    #
# REQUIRED:                                                                           #
#   - RefSeq exons                                                                    #
#   - RefSeq introns                                                                  #
#   - ${LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*             #
#                                                                                     # 
# 3) get_genesCoordinates.py                                                          #
# REQUIRED:                                                                           #
#   - ${LOVELACE}Annotation/Transcripts/grch37.clin.exons.refseq_ensembl_lrg_hugo_v*.bed.gz
#                                                                                     # 
#########################################################################################

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
##__________ SETUP __________## 
GENOMEDARCHIVE="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
if [ ! -d "${GENOMEDARCHIVE}" ]; then
  GENOMEDARCHIVE="/genomedarchive/"
fi
CRICK=${GENOMEDARCHIVE}Crick_storage/
LOVELACE=${GENOMEDARCHIVE}Lovelace_decoding/
MENDEL=${GENOMEDARCHIVE}Mendel_annotating/

############################################################
#   1) check for updates in grch37.clin.manual.refseq_ensembl.txt
############################################################

#if file doesn't exist, leave script
if [ ! -f ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update ]; then
    echo "File grch37.clin.manual.refseq_ensembl.txt doesn't have any updates."
    echo "There's no need to re-build RefSeq annotation files."
    exit
#if file exists, download it from GENEVA and update it
else
    echo "There's an update in grch37.clin.manual.refseq_ensembl.txt"
    echo "Getting grch37.clin.manual.refseq_ensembl.txt from GENEVA."
    mkdir ${MENDEL}clinical_manual/
    wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/annotations/grch37.clin.manual.refseq_ensembl.txt -P ${MENDEL}clinical_manual/

    #for each line in .update file, replace corresponding gene in grch37.clin.manual.refseq_ensembl.txt
    while IFS=$'\t' read -r -a line || [[ -n "$line" ]]; do
    echo "Replace or add gene ${line[0]} in grch37.clin.manual.refseq_ensembl.txt"
    grep -q "${line[0]}" ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.txt && sed -i "s/${line[0]}.*/${line[0]}\t${line[1]}\t${line[2]}\t${line[3]}\t${line[4]}/g" ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.txt || sed "$ a\ ${line[0]}\t${line[1]}\t${line[2]}\t${line[3]}\t${line[4]}" -i ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.txt
    done < ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update
fi

#remove .update file
rm ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update

#sort uniq grch37.clin.manual.refseq_ensembl.txt without header (body function)
less ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.txt | body sort -k1,1 -k2,2 - | uniq > ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.tmp
mv ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.tmp ${MENDEL}clinical_manual/grch37.clin.manual.refseq_ensembl.txt

############################################################
#   2) Run create_grch37.exons.refseq_ensembl_lrg_hugo.R
############################################################
${MENDEL}bin/create_grch37.refseq_ensembl_lrg_hugo.R

#______________________________________________
# version control
#______________________________________________
${MENDEL}bin/control_versions.sh

############################################################
#   3) Download BED files from UCSC Table Browser for create_grch37.refseq_ensembl_lrg_hugo.R
############################################################
# get exons.bed from UCSC Table Browser
wget --progress=dot \
'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&'\
'hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_compressType=none&outBed=1&'\
'fbQual=exon&fbExonBases=0&submit=submit&hgta_doGetBed=1' -O "${MENDEL}GRCh37_NCBIRefSeq_exons.bed"

# get introns.bed from UCSC Table Browser
wget --progress=dot \
'http://genome.ucsc.edu/cgi-bin/hgTables?clade=mammal&org=Human&db=hg19&'\
'hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=bed&hgta_compressType=none&outBed=1&'\
'fbQual=intron&fbIntronBases=0&submit=submit&hgta_doGetBed=1' -O "${MENDEL}GRCh37_NCBIRefSeq_introns.bed"

############################################################
#   4) Run create_grch37.exons.refseq_ensembl_lrg_hugo.R
############################################################
RefSeqtranscripts=$(ls ${LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*)

${MENDEL}/bin/create_grch37.exons.refseq_ensembl_lrg_hugo.R --exons=="${MENDEL}GRCh37_NCBIRefSeq_exons.bed" --introns=="${MENDEL}GRCh37_NCBIRefSeq_introns.bed" --RefSeq==${RefSeqtranscripts}
#______________________________________________
# version control
#______________________________________________
${MENDEL}bin/control_versions.sh

############################################################
#   4) Run get_genesCoordinates.py for clinical BED files
############################################################
RefSeqexons=$(ls ${LOVELACE}Annotation/Transcripts/grch37.clin.exons.refseq_ensembl_lrg_hugo_v*.bed.gz)

${MENDEL}/bin/get_genesCoordinates.py --refseq ${RefSeqexons} --out_file ${MENDEL}grch37.clin.exons.refseq_ensembl_lrg_hugo_coordinates.bed
#______________________________________________
# version control
#______________________________________________
${MENDEL}bin/control_versions.sh


