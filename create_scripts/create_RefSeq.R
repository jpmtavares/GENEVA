#!/usr/bin/Rscript

###########################################################################
#                                                                         #
#    Create NCBI RefSeq BED file (UCSC hg19) with clinical transcripts    #
#                                                                         #
###########################################################################
##!!! UCSC files are 0-based
#______________________________________________
# libraries
#______________________________________________
library(stringr)
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
if("--help" %in% args | is.null(args$exons) | is.null(args$introns) | is.null(args$clinical)) {
  cat("
      create_refSeq.R [--help] [--exons==<path to UCSC BED file with exons>] [--introns==<path to UCSC BED file with introns>] [--clinical==<path to file with clinical transcripts>]
      
      -- program to create BED file with RefSeq clinical transcripts along with ranking for exons and introns --
      
      where:
      --exons                      path to UCSC BED file with exons
      --introns                    path to UCSC BED file with introns
      --clinical                   path to file with clinical transcripts\n\n")
  
  q(save="no")
}

# SETUP
#________________________________________________________
exons<-read.delim(args[["exons"]], header=F)
introns<-read.delim(args[["introns"]], header=F)
clinical<-read.delim(args[["clinical"]], header = TRUE)
filename<-"refSeq_clinical"

################################################################
# 1) join exons and introns
################################################################
genes<-rbind(exons, introns)
names(genes)<-c("Chr", "Start", "End", "Info", "Zero", "Strand")

################################################################
# 2) split info
################################################################
genes<-data.frame(genes[,c("Chr", "Start", "End")],
                  str_split_fixed(genes$Info, "_", 8),
                  "Strand"=genes[,c("Strand")]) %>%
  mutate(Transcripts=paste(X1, X2, sep="_")) %>%
  select(Chr, Start, End, Transcripts, X3, X4, Strand)

names(genes)<-c("Chr", "Start", "End", "Transcripts", "Region", "Rank", "Strand")

################################################################
# 3) correct Region order in minus Strand and add clinical transcripts information
################################################################
bed<-genes %>%
  mutate(Start=Start+1) %>% #get 1-based document
  mutate(Region=ifelse(Region=="exon", "E", "I")) %>%
  group_by(Transcripts, Region) %>%
  mutate(Rank = ifelse(Strand=="-", # revert region order in minus Strand
                       rev(as.numeric(as.character(Rank))+1),
                       as.numeric(as.character(Rank))+1),
         Total=n(), # count total number of exons/introns in each transcript
         Chr=gsub("chr","", Chr)) %>% # remove "chr" 
  ungroup() %>%
  mutate(Region=paste(Region, Rank, "/", Total, sep="")) %>%
  mutate(refSeq_mRNA=str_extract(Transcripts, "[^\\.]*")) %>%
  inner_join(., clinical) %>%
  arrange(Chr, Start, End) %>%
  select(Chr, Start, End, Region, Strand, HGNC_symbol, ENSGene, ENSTranscript, refSeq_mRNA, refSeq_protein)

# turn-off scientific notation
bed$Start <-format(bed$Start, scientific = FALSE, trim=T) 
bed$End <-format(bed$End, scientific = FALSE, trim=T) 

################################################################
# 4) write-table
################################################################
write.table(bed, paste(filename, ".bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

################################################################
# 5) sort BED file, create BED for coverage analysis and indexed it
################################################################
system(paste("sort -k1,1 -k2,2n -k3,3n", paste(filename, ".bed", sep=""), ">", paste(filename, "_sort.bed", sep=""), sep=" "))
# BED for coverage analysis
system(paste("awk 'NF{NF-=1};1'", paste(filename, "_sort.bed", sep=""), "| uniq | awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"$6\",\"$4\",\"$5\",\"$7\",\"$8\",\"$9}' | sed -r 's/\\s+/\\t/g' >", paste(filename, "_coverage.bed", sep=""), sep=" "))
# index "_sort.bed"
system(paste("bgzip", paste(filename, "_sort.bed", sep=""), sep=" "))
system(paste("tabix -b 2 -e 3", paste(filename, "_sort.bed.gz", sep=""), sep=" "))
