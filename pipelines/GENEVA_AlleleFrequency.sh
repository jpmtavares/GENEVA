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

############################################################
#   1) Download get_scripts
############################################################
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_vcfmerge2freq.sh
wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/get_scripts/get_inHouseFreq.py

# Make them executables
chmod +x ${vcfmerge2freq}
chmod +x ${inHouseFreq}

############################################################
#   2) Run get_vcfmerge2freq.sh
############################################################
# Run with default parameters: present work directory will be used to get VCF files and write Allele frequency output file
./get_vcfmerge2freq.sh

############################################################
#   3) Remove get_scripts
############################################################
rm ${vcfmerge2freq}
rm ${inHouseFreq}
