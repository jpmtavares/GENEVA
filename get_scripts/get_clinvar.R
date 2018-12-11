#!/usr/local/bin/Rscript
#######################################################################################
#                                                                                     #
#                             ClinVar FTP download                                    #
#                                                                                     #
#######################################################################################
#' Get ClinVar GRCh37 vcf file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
#' 
#' @param d_ClinVar download directory A number.
#' @param y A number.
#' @return ClinVar GRCH37 vcf, remove older version, and replace date of updated file in config/vcfanno/myanno.conf 
#' @examples
#' Rscript ./get_clinvar.R
#' Rscript ./get_clinvar.R <ClinVar download directory> <vcfanno config file directory>
#______________________________________________
# setup: libraries
#______________________________________________
library(RCurl)
library(lubridate)
library(stringr)

#______________________________________________
# setup: variables
#______________________________________________
# passing arguments
args <- commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  args<-c("/usr/local/share/bcbio/genomes/Hsapiens/GRCh37/variation/",
          "/usr/local/share/bcbio/genomes/Hsapiens/GRCh37/config/vcfanno/")
} else if (length(args)!=2) {
  stop("Please insert 2 arguments: <download directory> and <vcfanno config file directory>.", call.=FALSE)
}

# get ClinVar vcf file
url<-"ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].vcf.gz" #doesn't accept [0-9]{8}
filename<-grep("clinvar_[0-9]*.vcf.gz$",unlist(strsplit(getURL(url, dirlistonly = TRUE), split = "\n")),
               value=T)

# get date of updated file in NCBI ClinVar
date_updated<-ymd(str_extract(filename,"[0-9]{8}"))
date_last<-ymd(max(str_extract(list.files(args[1], pattern = "clinvar*"),
                                 "[0-9]{8}")))
#______________________________________________
# download file
#______________________________________________
if(date_updated > date_last){ #if there is a newer file in NCBI Clinvar
    # then download it
    download.file(url, paste(args[1],filename,sep=""), method = "wget")
    # and index file
    download.file(paste(url, ".tbi", sep=""), paste(args[1],filename, ".tbi", sep=""), method = "wget")

    #remove older version
    file.remove(list.files(args[1], pattern = "clinvar*", full.names = TRUE)[1:2])

#______________________________________________
# update configuration file
#______________________________________________
conf<-readLines(paste(args[2], "myanno.conf", sep=""))
conf_updated<-gsub(pattern = "[0-9]{8}", replace = str_extract(filename,"[0-9]{8}"), x = conf) # replace older date by updated

writeLines(conf_updated, con=paste(args[2], "myanno.conf", sep=""))

}
