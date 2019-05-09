#!/usr/bin/env bash
#set -e 
#set -o pipefail

#######################################################################################
#                                                                                     #
#                               Control versions                                      #
#                                                                                     #
#######################################################################################


############################################################
#   HELP function
############################################################

usage="$(basename "$0") [-h]
       -- script that manages and controls file versions --
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
#   FUNCTIONS
############################################################

########## CREATES NEW ENTRY IF FILE DOESN'T EXIST ########## 
function create_entry {
    created_version="1.0"
    created_changes="created"
    echo
    echo "          CREATING FILE"
    echo
    echo "#####################################"
    echo "${filename_noVersion}"
    echo "#####################################"
    echo "version:" ${created_version}
    echo "date:" ${today}
    echo
    read -p "description [write a little description of the file]: "
    entry_description=$REPLY
    echo
    read_header #function that reads header of file and saves it to ${entry_header}
    echo "header:" ${entry_header}
    echo
    read -p "is_sorted [yes/no]: "
    entry_sort=$REPLY
    read -p "is_1based [yes/no]: "
    entry_1based=$REPLY
    echo
    echo "changes:" ${created_changes}
    echo

    echo -e "${filename_noVersion}\t${created_version}\t${today}\t${entry_description}\t"${entry_header}"\t${entry_sort}\t${entry_1based}\t${created_changes}" >> ${CRICK}Source_control/file_versions.info
  
    rename  "s/${filename_noVersion}/${filename_noVersion}_v${created_version}/" *
}

########## UPDATES ENTRY IF FILE ALREADY EXISTS ########## 
function add_entry {
    echo
    echo "          UPDATING FILE"
    echo
    echo "#####################################"
    echo "${filename_noVersion}"
    echo "#####################################"
    echo "version:" ${output_version}
    echo "date:" ${today}
    echo
    echo "description:" ${entry_fields[3]}
    echo
    echo "header:" ${entry_fields[4]}
    echo
    echo "is_sorted:" ${entry_fields[5]}
    echo "is_1based:" ${entry_fields[6]}
    echo
    read -p "changes [write what was changed]: "
    echo

    echo -e "${filename_noVersion}\t${output_version}\t${today}\t${entry_fields[3]}\t"${entry_fields[4]}"\t${entry_fields[5]}\t${entry_fields[6]}\t$REPLY" >> ${CRICK}Source_control/file_versions.info
  
    rename  "s/${filename_noVersion}/${filename_noVersion}_v${output_version}/" $(ls | grep "$filename_noVersion" | grep -v "create" | grep -v "get") 
}

########## SHOWS EXISTING ENTRY ##########
function existing_entry {
    echo
    echo "           EXISTING ENTRY"
    echo
    echo "#####################################"
    echo "${filename_noVersion}"
    echo "#####################################"
    echo "version:" ${entry_fields[1]}
    echo "date:" ${entry_fields[2]}
    echo
    echo "description:" ${entry_fields[3]}
    echo
    echo "header:" ${entry_fields[4]}
    echo
    echo "is_sorted:" ${entry_fields[5]}
    echo "is_1based:" ${entry_fields[6]}
    echo
    echo "changes:" ${entry_fields[7]}
    echo
    echo
}

########## READS HEADER OF ANNOTATION FILE ##########
function read_header {
    header=$(less $(ls | grep "^${filename_noVersion}" | grep -v "warnings") | head -n1)
    entry_header="${header//$'\t'/,}"
}

########## COPY FILES TO LOVELACE and CRICK ##########
function copy_files {
    crick_path=$1
    lovelace_path=$2

    # Copy new version to CRICK
    echo
    echo "Copy new file version to CRICK_storage/"
    echo
    cp $(ls | grep "^${filename_noVersion}_v" | grep -v "log\|warnings") ${crick_path}
    # Copy log and warnings files to CRICK
    cp $(ls | grep "${filename_noVersion}_v" | grep "log\|warning") ${crick_path}log/

    # Copy new version to LOVELACE and remove older one
    for link in $(ls | grep "^${filename_noVersion}_v" | grep -v "log\|warnings"); do
      echo
      echo "Remove older file version from LOVELACE_decoding/"
      echo
      rm ${lovelace_path}${link/_v[0-9]*.[0-9]/}
      echo
      echo "Copy new file version to LOVELACE_decoding/"
      echo
      ln -s ${crick_path}${link} ${lovelace_path}${link/_v[0-9]*.[0-9]/}
    done
    # Remove files in MENDEL
    rm $(ls | grep "${filename_noVersion}_v")
}

##__________ SETUP __________## 
GENOMEDARCHIVE="/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
if [ ! -d "${GENOMEDARCHIVE}" ]; then
  GENOMEDARCHIVE="/genomedarchive/"
fi
CRICK=${GENOMEDARCHIVE}Crick_storage/
LOVELACE=${GENOMEDARCHIVE}Lovelace_decoding/
##_______________________##
lognames=$(ls | grep "log_")
today=$(date '+%Y%m%d')
##_______________________##
for logfile in ${lognames}; do

filename="${logfile:4:-4}"
filenamelength=${#filename}
if [[ "${filename:$filenamelength-1:1}" == "." ]]; then
  filename="${filename::-1}"
fi

if [[ $filename =~ .*_v.* ]]; then
  filename_noVersion="${filename%_*}"
  version="${filename##*_v}"
else
  filename_noVersion=${filename}
fi
##_______________________##
entry=$(grep -w "${filename_noVersion}" ${CRICK}Source_control/file_versions.info | sort -k3n | tail -n1)
##_______________________##

############################################################
#   MAIN
############################################################
##########################
#   1) EVALUATE IF ENTRY EXISTS if file_versions.info
##########################
if [ -z "$entry" ]; then

    create_entry
    ##########################
    #   2) Copy to CRICK and LOVELACE
    ##########################
    read -p "Where in CRICK should I move new files?`echo $'\n> '`" -i "${CRICK}Annotation/" -e crick
    echo
    read -p "Where in LOVELACE should I move new files?`echo $'\n> '`" -i "${LOVELACE}Annotation/" -e lovelace
    echo
    copy_files ${crick} ${lovelace}

else
# create array entry_fields with ${entry} splitted by \t
IFS=$'\t' read -r -a entry_fields <<< "$entry"
  #entry_fields[0] filename
  #entry_fields[1] version
  #entry_fields[2] date
  #entry_fields[3] description
  #entry_fields[4] header
  #entry_fields[5] sort
  #entry_fields[6] 1based
  #entry_fiedls[7] changes
output_version=$(echo "${entry_fields[1]}+0.1" | bc)

  if [[ -z "${version}" ]]; then

    if [[ ${today} -eq ${entry_fields[2]} ]]; then
      echo
      echo "The last entry from this file is from today:"
      existing_entry
    
      read -p "Are you sure you want to update it? " -n 1 -r
      echo
    
      if [[ $REPLY =~ ^[Yy]$ ]];
      then
        add_entry
        ##########################
        #   2) Copy to CRICK and LOVELACE
        ##########################
        crick=$(dirname $(find ${CRICK} | grep ${filename_noVersion} | grep -v "log" | head -n1))/
        lovelace=$(dirname $(find ${LOVELACE} | grep ${filename_noVersion} | head -n1))/
        copy_files ${crick} ${lovelace}
      else
        exit 2
      fi
    
    else
      add_entry
      ##########################
      #   2) Copy to CRICK and LOVELACE
      ##########################
      crick=$(dirname $(find ${CRICK} | grep ${filename_noVersion} | grep -v "log" | head -n1))/
      lovelace=$(dirname $(find ${LOVELACE} | grep ${filename_noVersion} | head -n1))/
      copy_files ${crick} ${lovelace}
    fi
  else
    entry_existing=$(grep "$filename_noVersion" ${CRICK}Source_control/file_versions.info | grep -w ${version})
    if [[ -z ${entry_existing} ]]; then
      echo
      echo "!!!!!   Something went wrong   !!!!!"
      echo
      echo "Version ${version} of ${filename} doesn't exist in file_versions.info,"
      echo "please check this."
      echo
    else
      IFS=$'\t' read -r -a entry_fields <<< "$entry_existing"
        #entry_fields[0] filename
        #entry_fields[1] version
        #entry_fields[2] date
        #entry_fields[3] description
        #entry_fields[4] header
        #entry_fields[5] sort
        #entry_fields[6] 1based
        #entry_fields[7] changes

      existing_entry
      ##########################
      #   2) Copy to CRICK and LOVELACE
      ##########################
      crick=$(dirname $(find ${CRICK} | grep ${filename_noVersion} | grep -v "log" | head -n1))/
      lovelace=$(dirname $(find ${LOVELACE} | grep ${filename_noVersion} | head -n1))/
      copy_files ${crick} ${lovelace}

    fi
  fi
fi

done
