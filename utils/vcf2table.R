#!/usr/bin/Rscript

##########################################################
#                                                        #
#                      VCF2table                         #
#                                                        #
##########################################################
#______________________________________________
# libraries
#______________________________________________
library(bedr)
library(stringr)
library(magrittr)
options(scipen=999) # prevent scientific notation [to go back, set scipen to 0]
#______________________________________________
# HELP function
#______________________________________________
# collect arguments
args <- commandArgs(TRUE)

# parse arguments (in the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "==")
argsL <- as.list(as.character(as.data.frame(do.call("rbind", parseArgs(args)))$V2))

if(length(argsL)==0 & !is.null(argsL)){
  argsL <- as.list("--help")
}

names(argsL) <- as.data.frame(do.call("rbind", parseArgs(args)))$V1
args <- argsL
rm(argsL)

# Default setting when no all arguments passed or help needed
if("--help" %in% args | is.null(args$vcf)) {
  cat("
      vcf2table.R [--help] [--vcf==<path to VCF file>] [--header==<list of desired VCF fields>]

      -- converts a .VCF file into a .TSV --
      
      where:
      --vcf                        path to VCF file
      --header                     [default: all fields] list of desired VCF fields; last column of VCF must be called by 'SAMPLE'\n\n")
  
  q(save="no")
}

# SETUP
#________________________________________________________
sample<-str_extract(basename(args[["vcf"]]), ".*-..")
dir<-dirname(args[["vcf"]])

################################################################
# 1) read VCF file
################################################################
print("Reading table...")
myanno<-read.vcf(args[["vcf"]], split.info=T)
vcf<-myanno$vcf
# change last column name to SAMPLE
names(vcf)[ncol(vcf)]<-"SAMPLE"

################################################################
# 2) selection of the disered fields
################################################################
# if --header argument was not provided
if(length(args[["header"]])==0){
  header<-names(vcf) # header must contain all fields
}else{
  # if --header argument was provided
  header<-unlist(str_split(args[["header"]], "[, ]+|\t|\n")) %>% # the list is splitted consedering all delimeters: , \s \t \n 
    gsub("[, ()=]+|\n","",.) %>%
    .[. != ""] # blank elements are removed
}

# subset of selected fields
if(any(!(header %in% names(vcf)))){
  stop(paste("Selected field", header[!(header %in% names(vcf))], "is not in the VCF file.\nPlease remove it from list or change it by SAMPLE if you want the last VCF column.\n\n"))
} else {
  vcf<-vcf[header]
}

print(paste("Selected field: ", header, sep=" "))

################################################################
# 3) write.table
################################################################
print("Writing table...")
write.table(vcf, paste(dir, "/", sample,".table", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
print("Done!")
