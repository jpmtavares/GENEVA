#!/usr/bin/Rscript
#######################################################################################
#                                                                                     #
#                             ClinVar FTP download                                    #
#                                                                                     #
#######################################################################################
#' Get ClinVar GRCh37 vcf file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/
#' 
#' @param d_ClinVar download directory; default: /usr/local/share/bcbio/genomes/Hsapiens/GRCh37/variation/
#' @param d_vcfanno vcfanno config directory; default: /usr/local/share/bcbio/genomes/Hsapiens/GRCh37/config/vcfanno/
#' @return ClinVar GRCH37 vcf, remove older version, and replace date of updated file in config/vcfanno/myanno.conf 
#' @examples
#' Rscript ./get_clinvar.R
#' Rscript ./get_clinvar.R <ClinVar download directory> <vcfanno config file directory>
#______________________________________________
# libraries
#______________________________________________
library(RCurl, quietly = T)
library(lubridate, quietly = T, warn.conflicts = F)
library(stringr, quietly = T)
#__________________________________________________________
# SET GENOMEDARCHIVE and other PATHS
#__________________________________________________________
GENOMEDARCHIVE<-"/media/joanatavares/716533eb-f660-4a61-a679-ef610f66feed/"
if(!dir.exists(GENOMEDARCHIVE)){
  GENOMEDARCHIVE<-"/genomedarchive/"
}
CRICK<-paste(GENOMEDARCHIVE, "Crick_storage/", sep="")
LOVELACE<-paste(GENOMEDARCHIVE, "Lovelace_decoding/", sep="")
MENDEL<-paste(GENOMEDARCHIVE, "Mendel_annotating/", sep="")

setwd(MENDEL)
#______________________________________________
# HELP function
#______________________________________________
# collect arguments
args <- commandArgs(TRUE)

# parse arguments (in the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "==")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))

if(length(argsL)==0){
  args<-list()
  args$clinvarDownload<-paste(CRICK, "Annotation/Variants/ClinVar/", sep="")
  args$vcfannoConfig<-list.files(path=paste(LOVELACE, "config/", sep=""),
                                 pattern="grch37.vcfanno_germline_v*", full.names=T)
} else { 
  names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
  args <- argsL
  rm(argsL)
}

# Default setting when no all arguments passed or help needed
if("--help" %in% args) {
  cat("
      get_clinvar.R [--help] [--clinvarDownload==<dir to download clinvar>] [--vcfannoConfig==<path to vcfanno config file>]
      
      -- script to get ClinVar GRCh37 vcf file from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/ --
      
      where:
      --clinvarDownload            [default: ${CRICK}Annotation/Variants/ClinVar/] directory to download ClinVar\n\n")
  
  q(save="no")
}

#________________________________________________________
# SETUP
#________________________________________________________
#______________________________________
# START LOG for each file
#______________________________________
#sink(file = "log_grch37.vcfanno_germline.txt", append = FALSE, type = c("output"),
#     split = TRUE)

cat(paste("[", now("GMT"), "]", "Updating ClinVar", "\n"))

#______________________________________
# Read input files
#______________________________________
clinvar_d<-args[["clinvarDownload"]]
#vcfanno_d<-args[["vcfannoConfig"]]

#______________________________________________
# setup: variables
#______________________________________________
# passing arguments
#args <- commandArgs(trailingOnly=TRUE)

#if (length(args)==0) {
#  args<-c("/usr/local/share/bcbio/genomes/Hsapiens/GRCh37/variation/",
#          "/usr/local/share/bcbio/genomes/Hsapiens/GRCh37/config/vcfanno/")
#} else if (length(args)!=2) {
#  stop("Please insert 2 arguments: <download directory> and <vcfanno config file directory>.", call.=FALSE)
#}

# get ClinVar vcf file
url<-"ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar_[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].vcf.gz" #doesn't accept [0-9]{8}
filename<-grep("clinvar_[0-9]*.vcf.gz$",unlist(strsplit(getURL(url, dirlistonly = TRUE), split = "\n")),
               value=T)

# get date of updated file in NCBI ClinVar
date_updated<-ymd(str_extract(filename,"[0-9]{8}"))
date_last<-ymd(max(str_extract(list.files(clinvar_d, pattern = "clinvar*"),
                                 "[0-9]{8}")))

if(is.na(date_last)){
  date_last<-ymd(date_updated)-days(5)
}
#______________________________________________
# download file
#______________________________________________
if(date_updated > date_last){ #if there is a newer file in NCBI Clinvar
    # then download it
    cat(paste("                       ", "Downloading", filename, "\n"))
    download.file(url, paste(clinvar_d, filename, sep=""), method = "wget")
    # and index file
    download.file(paste(url, ".tbi", sep=""), paste(clinvar_d, filename, ".tbi", sep=""), method = "wget")

    #remove older versions from LOVELACE
    cat(paste("                       ", "Remove older version from LOVELACE", "\n"))
    file.remove(list.files(paste(LOVELACE, "Annotation/Variants/ClinVar/", sep =""), pattern = "clinvar*", full.names = TRUE)[1:2])
    #cp newer version to LOVELACE
    cat(paste("                       ", "Copy updated version to LOVELACE", "\n"))
    system(paste("cp", paste(clinvar_d, filename, sep=""), paste(LOVELACE, "Annotation/Variants/ClinVar/", sep ="")))
    system(paste("cp", paste(clinvar_d, filename, ".tbi", sep=""), paste(LOVELACE, "Annotation/Variants/ClinVar/", sep ="")))

    cat(paste("                       ", "ClinVar update finished with SUCCESS!!", "\n"))
#______________________________________________
# update configuration file
#______________________________________________
#system(paste("cp ", vcfanno_d, " ", MENDEL, "grch37.vcfanno_germline.conf", sep=""))
#conf<-readLines(paste(MENDEL, "grch37.vcfanno_germline.conf", sep=""))
#conf_updated<-gsub(pattern = "[0-9]{8}", replace = str_extract(filename,"[0-9]{8}"), x = conf) # replace older date by updated

#writeLines(conf_updated, con=paste(MENDEL, "grch37.vcfanno_germline.conf", sep=""))

} else {
  cat(paste("                       ", "The last version available on ClinVar FTP is from", date_updated, "which is the same as the last version found locally:", list.files(paste(LOVELACE, "Annotation/Variants", sep =""), pattern = "clinvar*", full.names = TRUE)[1], "\n"))
  cat(paste("                       ", "Therefore, there's no need to update ClinVar.", "\n"))
}
