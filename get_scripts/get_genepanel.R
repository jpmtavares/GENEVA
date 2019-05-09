#!/usr/bin/Rscript

###############################################################################
#                                                                             #
#          Get gene panel for each sample from fromana_yyyymmdd.tmp           #
#                                                                             #
###############################################################################
#                                                                             #
# This script reads a list of genes to be analyzed in each sample. It deals   #
# with some special characters and white-spaces that might exist in said list.#
#                                                                             #
# 1)                                                                          #
# This script reads a file with the format of                                 #
# ${CRICK}config/fromana_yyyymmdd.template that is stored in                  #
# ${CRICK}Source_control/gene_panel/ and outputs                              #
# ${samplename}_genepanel.txt in (1) ${CRICK}Raw/ || (2) ${LOVELACE}Raw/ ||   #
# (3) ${CRICK}Source_control/gene_panel/                                      #
#                                                                             #
# 2)                                                                          #
# If some genes in the list are not found in the REFERENCE ANNOTATION file    #
# ${LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*, it     #
# will output ${samplename}_genepanel.pending in the same directory.          #
#                                                                             #
# 3)                                                                          #
# If there is a gene in the list that doesn't have a clinical transcript,     #
# the user will be prompt about that and choose the appropriate one. It will  #
# then be outputed ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update  #
#                                                                             #
# NOTE: GENEVA_buildRefSeq.sh will evaluate if                                #
# ${LOVELACE}to_do/grch37.clin.manual.refseq_ensembl.update exists, and if    #
# so, it will re-build RefSeq annotations file with the clinical transcript   #
# update.                                                                     #  
#                                                                             #
# REQUIRED:                                                                   #
#   - fromana_yyyymmdd.tmp [--fromana]                                        #
#   - ${LOVELACE}Annotation/Transcripts/grch37.refseq_ensembl_lrg_hugo_v*     # 
#                                                                             #
###############################################################################


#______________________________________________
# libraries
#______________________________________________
#install.packages("dplyr")
#install.packages("lubridate")
library(dplyr, quietly = T, warn.conflicts = F)
library(lubridate, quietly = T, warn.conflicts = F)
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

#__________________________________________________________
# HELP function
#__________________________________________________________
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
if("--help" %in% args | is.null(args$fromana)) {
  cat("
      get_genepanel.R [--help] [--fromana==<path to file sent by AC>]
      -- program that gets a list of genes that will be analysed in each sample, and confirms
         if the list is well written and has clinical transcripts associated --
      
      where:
      --fromana       [default: Crick_storage/Source_control/gene_panel/] path to file sent by AC fromana_yyyymmdd.tmp\n\n")
  
  q(save="no")
}
#__________________________________________________________
# SETUP
#__________________________________________________________
# default values
#________________________
if(any(grep("/", args$fromana))){
  fromana<-read.delim(args[["fromana"]], header=F)
} else {
  fromana<-read.delim(paste(CRICK, "Source_control/gene_panel/",args$fromana, sep=""),
                      header=F)
}

#__________________________________________________________
# FUNCTIONS
#__________________________________________________________

########## GENES OF INTEREST ########## 
genesofinterest<-function(sampleID, genesList){
  #________________________________________________________
  # read RefSeq transcripts
  #________________________________________________________
  refseq_path<-list.files(path=paste(LOVELACE, "Annotation/Transcripts/", sep=""),
                          pattern="grch37.refseq_ensembl_lrg_hugo_v*", recursive = F, full.names = T)
  refseq<-read.delim(refseq_path, header=T)
  
  #________________________________________________________
  # remove white spaces or other characters from genesList and sampleID
  #________________________________________________________
  # remove non-ASCII characters
  genesList<-iconv(genesList, "latin1", "ASCII", sub="")
  #processing genes list
  #_______________________
  genes<-genesList %>%
    gsub("\\(.*\\)","", .) %>% #removes everything between ( )
    gsub("\\[.*\\]","", .) #removes everything between [ ]
  genes<-unlist(strsplit(genes,split="[, \t=\\*]+")) %>% #splits list by ,\s + \t or =
    gsub("[, ()=]+|\n","",.) %>%
    .[.!=""] %>% unique %>% sort %>% data.frame(HGNC_symbol=.)
  #_______________________
  #processing sample IDs
  #_______________________
  sample<-gsub("[, ()=]+|\n","",sampleID)
  
  #______________________________________
  # START LOG for each sample
  #______________________________________
  sink(file = paste(LOVELACE, "log/log_",sample, sep=""),
       append = TRUE, type = c("output"), split = FALSE)
  #_______________________
  # write STEP and count genes for each sample
  #_______________________
  cat(paste("[", now("GMT"), "]", "Starting Gene Panel for", sample,"\n\n"))
  cat(paste("                       ", "There are", nrow(genes), "gene(s) in", basename(args$fromana),":", paste(genes$HGNC_symbol, collapse=", "), "\n\n", sep=" "))
  
  #________________________________________________________
  # check if genesList is well-written and has clinical transcripts
  #________________________________________________________
  hgnc_symbol<-left_join(genes, refseq) %>%
    group_by(HGNC_symbol) %>%
    filter(if(any(grep("yes", clinical_transcript))) clinical_transcript!="no" else clinical_transcript=="no")
  
  hgnc_alternative_symbol<-left_join(genes, refseq) %>%
    filter(is.na(HGNC_alternative_symbol)) %>%
    transmute(HGNC_alternative_symbol=HGNC_symbol) %>%
    left_join(., refseq) %>%
    group_by(HGNC_symbol) %>%
    filter(if(any(grep("yes", clinical_transcript))) clinical_transcript!="no" else clinical_transcript=="no") %>%
    select(HGNC_symbol, HGNC_alternative_symbol, refSeq_mRNA, refSeq_protein, refSeq_mRNA_noVersion, refSeq_protein_noVersion,
           ENSGene, ENSTranscript, LRG_id, clinical_transcript)
  
  bad_written<-left_join(genes, refseq) %>%
    filter(is.na(HGNC_alternative_symbol)) %>%
    transmute(HGNC_alternative_symbol=HGNC_symbol) %>%
    left_join(., refseq) %>%
    filter(is.na(HGNC_symbol)) %>%
    transmute(HGNC_symbol=HGNC_alternative_symbol)
  
  #________________________________________________________
  # genes with no clinical transcript
  #________________________________________________________
  if(any(hgnc_symbol$clinical_transcript=="no")){
    hgnc_symbol<-chooseclinical(hgnc_symbol, LOVELACE, sample)
  }
  if(any(hgnc_alternative_symbol$clinical_transcript=="no")){
    hgnc_alternative_symbol<-chooseclinical(hgnc_alternative_symbol, LOVELACE, sample)
  }
  
  full.info<-rbind(hgnc_symbol, hgnc_alternative_symbol) %>%
    unique %>%
    filter(clinical_transcript!="no")
  #________________________________________________________
  # print accordingly
  #________________________________________________________
  total_genes<-nrow(genes)
  #print all manually curated transcripts
  print_manual(full.info, total_genes)
  #IF not all transcripts were manually curated
  #they might be written in
  ############################### ALTERNATIVE FORM  ############################### 
  if(nrow(hgnc_alternative_symbol)>0){ 
    print_alternative(hgnc_alternative_symbol, total_genes)
    
    #they might have clinical transcript because they were chosen from user
    if(any(hgnc_alternative_symbol$clinical_transcript=="user_choice")){
      print_choice(hgnc_alternative_symbol, total_genes)
    }
    #they might have clinical transcript because there is only 1 transcript
    if(any(hgnc_alternative_symbol$clinical_transcript=="yes_automatic")){ 
      print_automatic(hgnc_alternative_symbol, total_genes)
    }
    #they might have clinical transcript because there is only 1 coding transcript
    if(any(hgnc_alternative_symbol$clinical_transcript=="yes_coding")){
      print_coding(hgnc_alternative_symbol, total_genes)
    }
    #they might have clinical transcript because they were obtained from HUGO genes
    if(any(hgnc_alternative_symbol$clinical_transcript=="yes_HUGO")){
      print_HUGO(hgnc_alternative_symbol, total_genes)
    }
  }
  ################################################################################# 
  #they might have clinical transcript because they were chosen from user
  if(any(hgnc_symbol$clinical_transcript=="user_choice")){
    print_choice(hgnc_symbol, total_genes)
  }
  #they might have clinical transcript because there is only 1 transcript
  if(any(hgnc_symbol$clinical_transcript=="yes_automatic")){ 
    print_automatic(hgnc_symbol, total_genes)
  }
  #they might have clinical transcript because there is only 1 coding transcript
  if(any(hgnc_symbol$clinical_transcript=="yes_coding")){
    print_coding(hgnc_symbol, total_genes)
  }
  #they might have clinical transcript because they were obtained from HUGO genes
  if(any(hgnc_symbol$clinical_transcript=="yes_HUGO")){
    print_HUGO(hgnc_symbol, total_genes)
  }
  #they might be bad-written or not included in RefSeq annotation
  if(nrow(bad_written)>0){ 
    print_badgene(sample, genesList, total_genes, bad_written, refseq_path, CRICK)
  }
  #______________________________________
  # write.file
  #______________________________________
  cat(paste("\n", "                        ", "OUTPUT GENE PANEL: ",
            paste(unique(sort(c(hgnc_symbol$HGNC_symbol, hgnc_alternative_symbol$HGNC_alternative_symbol))), collapse = ", "),
            "\n",sep=""))
  genes_output<-data.frame(genes=sort(c(hgnc_symbol$HGNC_symbol, hgnc_alternative_symbol$HGNC_alternative_symbol))) %>%
    unique
  
  # write file in CRICK (means that all went well in pre-processing)
  tryCatch(expr = {
    write.table(genes_output, paste(CRICK, "Raw/2019/",sample, "/", sample, "_genepanel.txt", sep = ""),
                col.names = F, row.names = F, quote = F)
  },
  # try to write file in LOVELACE (means that something went wrong during pre-processing)
  error = function(e) {
    write.table(genes_output, paste(LOVELACE, "Raw/",sample, "/", sample, "_genepanel.txt", sep = ""),
                col.names = F, row.names = F, quote = F)
  },
  # leave file in CRICK/Source_control/gene_panel (means that didn't find sample in LOVELACE)
  error = function(e2) {
    write.table(genes_output, paste(CRICK, "Source_control/gene_panel/", sample, "_genepanel.txt", sep = ""),
                col.names = F, row.names = F, quote = F)
    cat(paste("[        ERROR!!      ]", "Gene Panel failed for", "\n", sample, sep = " "))
    cat(paste("                       ", "This sample wasn't found in CRICK/Raw or LOVELACE/Raw... please check this.", "\n"))
  })
  
  sink()
}

########## CHOOSE CLINICAL TRANSCRIPTS ########## 
chooseclinical<-function(HGNC_genes, LOVELACE, sample){
  #______________________________________
  # STOP LOG to show question to user
  #______________________________________
  sink()
  #______________________________________
  # Read genes with no clinical transcript
  #______________________________________
  noclinical<-data.frame(HGNC_genes[grep("no", HGNC_genes$clinical_transcript),])
  cat("\n\n________________________________________________________________________________________\n\n")
  cat(paste("Gene", unique(noclinical$HGNC_symbol), "doesn't have a clinical transcript attributed.\n"))
  cat("Which of the following transcripts do you want to use?\n")
  cat("________________________________________________________________________________________\n\n")
  print(apply(noclinical,2,paste,sep="\t"))
  cat("________________________________________________________________________________________\n")
  cat(paste("Choose one number ", "[", 1, "-", nrow(noclinical),"] :", "\n", sep=""))
  answer<-readLines("stdin", n=1);
  #answer<-readline()
  #______________________________________
  # write answer to .update file
  #______________________________________
  write.table(noclinical[as.numeric(answer),c("HGNC_symbol", "ENSGene", "ENSTranscript",
                                              "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")],
              paste(LOVELACE, "to_do/", "grch37.clin.manual.refseq_ensembl.update", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
  
  #______________________________________
  # RESTART LOG
  #______________________________________
  sink(file = paste(LOVELACE, "log/log_",sample, sep=""),
       append = TRUE, type = c("output"), split = FALSE)
  
  HGNC_genes$clinical_transcript<-factor(HGNC_genes$clinical_transcript, levels=c(levels(HGNC_genes$clinical_transcript), "user_choice"))
  HGNC_genes$clinical_transcript[HGNC_genes$HGNC_symbol==noclinical$HGNC_symbol[as.numeric(answer)]][as.numeric(answer)] <- "user_choice"
  #HGNC_genes$clinical_transcript[as.numeric(answer)] <- "user_choice"
  
  return(HGNC_genes)
}

########## PRINT FUNCTIONS ##########
print_manual<-function(HGNC_genes, total_genes){
  manual<-unique(HGNC_genes$HGNC_symbol[HGNC_genes$clinical_transcript %in% "yes_manual"])
  cat(paste("                       ", "MANUALLY CURATED",
            paste("[", length(table(manual)), "/", total_genes,"]", sep=""),
            paste(manual, collapse=", "), "\n", sep=" "))
}
print_choice<-function(HGNC_genes, total_genes){
  choice<-unique(HGNC_genes$HGNC_symbol[HGNC_genes$clinical_transcript %in% "user_choice"])
  transcript<-unique(HGNC_genes$refSeq_mRNA[HGNC_genes$clinical_transcript %in% "user_choice"])
  cat(paste("                       ", "USER CHOICE",
            paste("[", length(table(choice)), "/", total_genes,"]", sep=""),
            paste(choice, collapse=", "),
            paste("(", paste(transcript, collapse=", "), ")",sep=""), "\n", sep=" "))
}
print_automatic<-function(HGNC_genes, total_genes){
  automatic<-unique(HGNC_genes$HGNC_symbol[HGNC_genes$clinical_transcript %in% "yes_automatic"])
  cat(paste("                       ", "AUTOMATIC-1 TRANSCRIPT",
            paste("[", length(table(automatic)), "/", total_genes,"]", sep=""),
            paste(automatic, collapse=", "),"\n", sep=" "))
}
print_coding<-function(HGNC_genes, total_genes){
  coding<-unique(HGNC_genes$HGNC_symbol[HGNC_genes$clinical_transcript %in% "yes_coding"])
  cat(paste("                       ", "AUTOMATIC-1 CODING",
            paste("[", length(table(coding)), "/", total_genes,"]", sep=""),
            paste(coding, collapse=", "),"\n", sep=" "))
}
print_HUGO<-function(HGNC_genes, total_genes){
  hugo<-unique(HGNC_genes$HGNC_symbol[HGNC_genes$clinical_transcript %in% "yes_HUGO"])
  cat(paste("                       ", "HUGO GENES",
            paste("[", length(table(hugo)), "/", total_genes,"]", sep=""),
            paste(hugo, collapse=", "),"\n", sep=" "))
}
print_badgene<-function(sample, genesList, total_genes, bad_written, refseq_path, CRICK){
  cat(paste("[        ERROR!!      ]", "Gene Panel failed for", sample, sep = " "))
  cat(paste("                       ", "GENES NOT INCLUDED in", basename(refseq_path),
            paste("[", length(bad_written$HGNC_symbol), "/", total_genes,"]", sep=""),
            paste(bad_written$HGNC_symbol, collapse = ", "), "\n", sep=" "))
  cat(paste("                       ", "PLEASE CHECK! Is there any typo in gene(s) name(s)?\n"))
  cat(paste("                       ", "File", paste(sample, "genepanel.pending", sep = "_"),"will be created.\n"))
  
  # write .pending file in CRICK (means that all went well in pre-processing)
  tryCatch(expr = {
    write.table(data.frame(sample, genesList),paste(CRICK, "Raw/2019/", sample, "/", sample, "_genepanel.pending", sep = ""),
                col.names = F, row.names = F, quote = F, sep="\t")
  },
  # try to write .pending file in LOVELACE (means that something went wrong during pre-processing)
  error = function(e) {
    write.table(data.frame(sample, genesList),paste(LOVELACE, "Raw/", sample, "/", sample, "_genepanel.pending", sep = ""),
                col.names = F, row.names = F, quote = F, sep="\t")
  },
  # leave .pending file in CRICK/Source_control/gene_panel (means that didn't find sample in LOVELACE)
  error = function(e2) {
    write.table(data.frame(sample, genesList),paste(CRICK, "Source_control/gene_panel/", sample, "_genepanel.pending", sep = ""),
                col.names = F, row.names = F, quote = F, sep="\t")
    cat(paste("[        ERROR!!      ]", "Gene Panel failed for", sample, sep = " "))
    cat(paste("                       ", "This sample wasn't found in CRICK/Raw or LOVELACE/Raw... please check this."))
  })
  
  # additionally, write .pending file in LOVELACE/to_do
  write.table(data.frame(sample, genesList),paste(LOVELACE, "to_do/", sample, "_genepanel.pending", sep = ""),
              col.names = F, row.names = F, quote = F, sep="\t")
    #stop 
  stop("ERROR!")
}
print_alternative<-function(HGNC_genes, total_genes){
  cat(paste("                       ", "ALTERNATIVE NAME                   |",
            paste("[", length(unique(HGNC_genes$HGNC_alternative_symbol)), "/", total_genes,"]", sep=""),
            paste(HGNC_genes$HGNC_alternative_symbol, collapse = ", "),
            paste("(", paste(HGNC_genes$HGNC_symbol, collapse = ", "), ")", sep=""),"\n", sep=" "))
}
#________________________________________________________
# MAIN
#________________________________________________________
#########################################################
# 1) Get samples IDs and Genes
#########################################################
fromana_sample_genes<-fromana[-c(1,2),] # removes download link and platform (2 first lines)
# apply genesofinterest to each line in fromana_yyyymmdd.tmp
apply(fromana_sample_genes, 1, function(x) {
  tryCatch({
    genesofinterest(as.character(x[1]),
                    as.character(x[2]))
  }, error = function(e){})  #if some gene name is not found
})                           #this will avoid the script to crash
#and will continue for the next sample

