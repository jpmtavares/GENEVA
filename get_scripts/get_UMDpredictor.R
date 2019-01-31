#!/usr/bin/Rscript

##########################################################
#                                                        #
#                     UMD-predictor                      #
#                                                        #
##########################################################
#______________________________________________
# libraries
#______________________________________________
library(readr)
library(stringr)
library(biomaRt)
library(dplyr)
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
if("--help" %in% args | is.null(args$ENSTranscripts)) {
  cat("
      get_UMDpredictor.R [--help] [--ENSTranscripts==<list of Ensembl transcripts ID>]

      -- program to get umd-predictor.eu prediction scores for a list of Ensembl Transcripts IDs --
      
      where:
      --ENSTranscripts                    list of Ensembl transcripts IDs\n\n")
  
  q(save="no")
}

#__________________________________________________________
# FUNCTIONS
#__________________________________________________________
UMDpredictor<-function(ENSTranscripts){
  # create url for each transcript
  url<-paste("http://umd-predictor.eu/transcript_query2.php?name=",
             str_extract(ENSTranscripts,"(?<=0{5}).*$"),sep="")
  # read trnascript information
  umd<-lapply(url, function(x) {
    tryCatch({
      print(paste("Get ENST00000", str_extract(x,"\\d{6}"), " transcript.", sep=""))
      read.delim(x, header = F)},
      error=function(e){NULL})
  })
  
  # remove error entries
  print(paste("Couldn't get results for", ENSTranscripts[sapply(umd, is.null)], "transcript.", sep=" "))
  ENSTranscripts<-ENSTranscripts[!sapply(umd, is.null)]
  umd<-umd[!sapply(umd, is.null)]
  
  # add ENSTranscript column
  umd<-mapply(cbind, "ENSTranscript"=ENSTranscripts, umd,  SIMPLIFY=F)
  # list 2 data.frame
  umdpredictor<-do.call("rbind", umd)
  names(umdpredictor)[-1]<-c("TranscriptPosition","Position","HGVS_c","HGVSp",
                             "score","UMD-predictor")
  #________________________________
  # Get info from Ensembl biomart
  #________________________________
  print("Getting Ensembl transcript information from BioMart...")
  ensembl<-useMart("ensembl", dataset="hsapiens_gene_ensembl")
  Ensembl_info<-getBM(attributes = c('chromosome_name', 'hgnc_symbol', 'ensembl_transcript_id'), 
                      filters = 'ensembl_transcript_id', 
                      values = ENSTranscripts, 
                      mart = ensembl)
  names(Ensembl_info)<-c("Chr", "HGNC_symbol", "ENSTranscript")
  print("BioMart retrieved.")
  
  # get REF and ALT alleles
  alleles<-do.call(rbind.data.frame, str_extract_all(as.character(umdpredictor$HGVS_c), "[ACGT]"))
  # join Ensembl and UMD-predictor
  print("Joining results...")
  result<-left_join(umdpredictor, Ensembl_info)
  result<-data.frame(Chr=result$Chr, Position=result$Position, Ref=alleles[,1], Alt=alleles[,2],
                     HGVS_c=result$HGVS_c, HGVS_p=result$HGVSp, HGNC_symbol=result$HGNC_symbol,
                     ENSTranscript=result$ENSTranscript, UMD_pred=result$`UMD-predictor`,
                     UMD_score=result$score)
  
  return(result)
}

#________________________________________________________
# start log file
#________________________________________________________
sink(file = "log.txt", append = FALSE, type = c("output"),
     split = FALSE)
#________________________________________________________
# SETUP
#________________________________________________________
transcripts<-read.delim(args[["ENSTranscripts"]], header=F)
filename<-"UMD-predictor_clinical_transcripts"

#########################################################
# 1) Check for strange Ensembl IDs
#########################################################
# if there is some row that doesn't contain an Ensembl ID
if(any(!grepl("ENST\\d{11}$", transcripts[,1]))) { 
  print(paste("Ignoring line", which(!grepl("ENST\\d{11}$", transcripts[,1])),
              "-", transcripts[which(!grepl("ENST\\d{11}$", transcripts[,1])),], "-",
              "from input file.", sep=" "))
  
  # Remove strange IDs from data.frame and get unique entries
  # BUT!! they are kept in input file
  transcripts<-as.character(unique(transcripts[-which(!grepl("ENST\\d{11}$", transcripts[,1])),]))
}

print(paste("Reading", length(transcripts), "transcripts.", sep=" "))

#########################################################
# 2) run function with corrected Ensembl IDs
#########################################################
umd<-UMDpredictor(transcripts)

#########################################################
# 3) write-table
#########################################################
print("Writing final table to file...")
write.table(umd, paste(filename,".txt",sep=""), col.names=F, row.names=F, quote=F, sep="\t")
sink() # close log file

#########################################################
# 4) sort otuput file and indexed it
#########################################################
system(paste("sort -k1,1 -k2,2n ", paste(filename,".txt",sep=""), ">", paste(filename, "_sort.txt", sep=""), sep=" "))
# index "_sort.txt"
system(paste("bgzip", paste(filename, "_sort.txt", sep=""), sep=" "))
system(paste("tabix -b 2 -e 2", paste(filename, "_sort.txt.gz", sep=""), sep=" "))
