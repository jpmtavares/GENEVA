#!/usr/bin/Rscript

##########################################################
#                                                        #
#                       Read VCF                         #
#                                                        #
##########################################################
#______________________________________________
# libraries
#______________________________________________
library(bedr)
library(stringr)
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
      vcf2table.R [--help] [--vcf==<path to VCF file>]

      -- converts a .VCF file into a .TSV --
      
      where:
      --vcf                        path to VCF file\n\n")
  
  q(save="no")
}

# SETUP
#________________________________________________________
print("Reading table...")
myanno<-read.vcf(args[["vcf"]], split.info=T)
sample<-str_extract(basename(args[["vcf"]]), ".*-..")
dir<-dirname(args[["vcf"]])

################################################################
# 1) write.table
################################################################
vcf<-myanno$vcf
print("Writing table...")
write.table(vcf, paste(dir, "/", sample,".table", sep=""), col.names = T, row.names = F, quote=F, sep="\t")
print("Done!")
