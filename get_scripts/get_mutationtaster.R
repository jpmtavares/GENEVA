#!/usr/bin/Rscript

##########################################################
#                                                        #
#                     MutationTaster                     #
#                                                        #
##########################################################

library(stringr, quietly=TRUE)
library(doParallel, quietly=TRUE)
library(XML, quietly=TRUE)
library(rlist, quietly=TRUE)
library(reticulate, quietly=TRUE)
use_python("/usr/bin/python3.5")

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
if("--help" %in% args | is.null(args$variants)) {
  cat("
      get_mutationtaster.R [--help] [--variants==<path to variants file with Chr, Position, Ref, Alt>]

      -- program to get mutationtaster.org scores --
      
      where:
      --variants                    path to variants file\n\n")
  
  q(save="no")
}

# SETUP
#________________________________________________________
variants<-read.delim(args[["variants"]], header=F)
filename<-basename(args[["variants"]])

names(variants)<-c("Chr", "Position", "Ref", "Alt")

#_________________________________________________________
# FUNCTIONS
#_________________________________________________________
################################################################
#              get mutationTaster - main function              # 
################################################################
MT<-function(variants){
  ##__________________________________
  ## parallelize
  ##__________________________________
  registerDoParallel(cores=1)
  print(paste("Getting MutationTaster scores for ", filename, "..."))
  parallel<-foreach(v=variants$Chr, w=variants$Position, x=variants$Ref,
                    y=variants$Alt) %dopar% data.frame(Chr=v, Position=w, Ref=x, Alt=y, mutationTaster(v,w,x,y))
  
  mt<-do.call(rbind.data.frame,parallel)
  print(paste("MutationTaster scores for ", filename, " retrieved."))

  write.table(mt[ ,c(1:4, 14, 5:7, 10)], paste("MT_", filename, sep="") , col.names=T, row.names=F, quote=F, sep="\t")
  system(paste("mv ", filename, " ./MT_success_files/", sep=""))
}

################################################################
#           Make a GET request to get mutationTaster           #
################################################################
mutationTaster<-function(Chr,Position,Ref,Alt){
  requests<-import("requests")
  parameters<-dict(chromosome=as.character(Chr),
                   position=as.character(Position),
                   ref=as.character(Ref),
                   alt=as.character(Alt))

  response<-requests$get("http://www.mutationtaster.org/cgi-bin/MutationTaster/MT_ChrPos.cgi?", params = parameters)

  #if there is some error with the query
  if(grepl("InDels must start with the last reference base|wrong input format|data problem|no suitable transcript|annotation problem|Software error",
            response$content)){
    problem<-str_extract(response$content,"InDels must start with the last reference base|wrong input format|data problem|no suitable transcript")
    results<-data.frame(genesymbol=NA, prediction=NA, probability=NA,
                        model=NA, predictionproblem=NA, splicing=NA,
                        ClinVar=NA, `amino acid changes`=NA,
                        `variant type`=NA, `dbSNP ID`=NA, `protein length`=NA,
                        file=NA)
  } else {
      tables<-readHTMLTable(response$content)
      tables <- list.clean(tables, fun = is.null, recursive = FALSE)
      results<-as.data.frame(tables[[1]])
  }
  
  return(unique(results))
}

#_________________________________________________________
# MAIN
#_________________________________________________________
sink(file = "log.txt", append = TRUE, type = c("output"),
     split = FALSE)
dir.create("./MT_success_files/", showWarnings = FALSE, recursive = FALSE)
print(paste("Reading ", filename, "...", sep=""))
# Run script
MT(variants)
print(paste(filename, "finished with SUCCESS!!", sep=" "))
sink()
