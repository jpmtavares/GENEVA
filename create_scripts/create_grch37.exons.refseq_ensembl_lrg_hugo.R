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
library(stringr, quietly = T)
library(plyr, quietly = T)
library(dplyr, quietly = T, warn.conflicts = F)
library(magrittr, quietly = T)
library(ape, quietly = T)
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

setwd(MENDEL)
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
if("--help" %in% args | is.null(args$exons) | is.null(args$introns) | is.null(args$RefSeq)) {
  cat("
      create_grch37.exons.refseq_ensembl_lrg_hugo.R [--help] [--exons==<path to UCSC BED file with exons>] [--introns==<path to UCSC BED file with introns>] [--RefSeq==<path to file with custom RefSeq annotation>]
      
      -- program to create BED file with RefSeq clinical transcripts along with ranking for exons and introns --
      
      where:
      --exons                      path to UCSC BED file with exons
      --introns                    path to UCSC BED file with introns
      --RefSeq                     path to file with custom RefSeq annotation\n\n")
  
  q(save="no")
}
#________________________________________________________
# SETUP
#________________________________________________________
#______________________________________
# START LOG for each file
#______________________________________
sink(file = "log_grch37.exons.refseq_ensembl_lrg_hugo.txt", append = FALSE, type = c("output"),
     split = TRUE)
sink(file = "log_grch37.clin.exons.refseq_ensembl_lrg_hugo.txt", append = FALSE, type = c("output"),
     split = TRUE)
sink(file = "log_grch37.clin.exons.refseq_ensembl_lrg_hugo_coverage.txt", append = FALSE, type = c("output"),
     split = TRUE)

cat("########################################################\n")
cat("create_grch37.exons.refseq_ensembl_lrg_hugo.R\n")
cat("########################################################\n")
cat(paste("[", now("GMT"), "]", "Creating GRCh37 exons/introns RefSeq annotation files", "\n"))

#______________________________________
# Read input files
#______________________________________
exons<-read.delim(args[["exons"]], header= F)
introns<-read.delim(args[["introns"]], header= F)
RefSeq<-read.delim(args[["RefSeq"]], header = TRUE)
filename<-"grch37.exons.refseq_ensembl_lrg_hugo"
filename_clinical<-"grch37.clin.exons.refseq_ensembl_lrg_hugo"

#______________________________________
# print input files
#______________________________________
cat(paste("                        ", "Input files: ", basename(args[["exons"]]),"\n",
"                                     ", basename(args[["introns"]]), "\n",
"                                     ", basename(args[["RefSeq"]]), "\n", sep=""))

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
  dplyr::select(Chr, Start, End, Transcripts, X3, X4, Strand)

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
  mutate(refSeq_mRNA_noVersion=str_extract(Transcripts, "[^\\.]*")) %>%
  inner_join(., RefSeq) %>%
  arrange(Chr, Start, End) %>%
  dplyr::select(Chr, Start, End, Region, Strand, HGNC_symbol, HGNC_alternative_symbol, ENSGene, ENSTranscript, refSeq_mRNA, refSeq_protein,
         refSeq_mRNA_noVersion, refSeq_protein_noVersion, LRG_id, clinical_transcript)

# get only clinical transcripts
bed_clinical<-bed %>%
  filter(clinical_transcript!="no") %>%
  dplyr::select(-clinical_transcript)

# turn-off scientific notation
bed$Start <-format(bed$Start, scientific = FALSE, trim=T) 
bed$End <-format(bed$End, scientific = FALSE, trim=T) 

# turn-off scientific notation
bed_clinical$Start <-format(bed_clinical$Start, scientific = FALSE, trim=T) 
bed_clinical$End <-format(bed_clinical$End, scientific = FALSE, trim=T) 
################################################################
# 4) write-table
################################################################
write.table(bed, paste(filename, ".bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
write.table(bed_clinical, paste(filename_clinical, ".bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

################################################################
# 5) sort BED file, create BED for coverage analysis and indexed it
################################################################
system(paste("sort -k1,1 -k2,2n", paste(filename, ".bed", sep=""), ">", paste(filename, "_sort.bed", sep=""), sep=" "))
system(paste("sort -k1,1 -k2,2n", paste(filename_clinical, ".bed", sep=""), ">", paste(filename_clinical, "_sort.bed", sep=""), sep=" "))

# BED for coverage analysis: Chr    Start    End    Gene,Region,Strand,ENSGene,ENSTranscript,RefSeqNM
#system(paste("awk 'NF{NF-=1};1'", paste(filename, "_sort.bed", sep=""), "| uniq | awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"$6\",\"$7\",\"$4\",\"$5\",\"$8\",\"$9\",\"$10}' | sed -r 's/\\s+/\\t/g' >", paste(filename, "_coverage.bed", sep=""), sep=" "))
system(paste("awk 'NF{NF-=1};1'", paste(filename_clinical, "_sort.bed", sep=""), "| uniq | awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"$6\",\"$7\",\"$4\",\"$5\",\"$8\",\"$9\",\"$10}' | sed -r 's/\\s+/\\t/g' >", paste(filename_clinical, "_coverage.bed", sep=""), sep=" "))

# add header to output files
write.table("#Chr\tStart\tEnd\tRegion\tStrand\tHGNC_symbol\tHGNC_alternative_symbol\tENSGene\tENSTranscript\trefSeq_mRNA\trefSeq_protein\trefSeq_mRNA_noVersion\trefSeq_protein_noVersion\tLRG_id", "header.txt", col.names=F, row.names=F, quote=F, sep="\t")
system(paste("cat header.txt", paste(filename, "_sort.bed", sep=""), ">",
             paste(filename, "_hdr.bed", sep=""), sep=" "))
system(paste("rm", paste(filename, "_sort.bed", sep="")))
system(paste("cat header.txt", paste(filename_clinical, "_sort.bed", sep=""), ">",
             paste(filename_clinical, "_hdr.bed", sep=""), sep=" "))
system(paste("rm", paste(filename_clinical, "_sort.bed", sep="")))

# bgzip ".bed"
system(paste("mv", paste(filename, "_hdr.bed", sep=""), paste(filename, ".bed", sep="")))
system(paste("bgzip", paste(filename, ".bed", sep=""), sep=" "))
system(paste("mv", paste(filename_clinical, "_hdr.bed", sep=""), paste(filename_clinical, ".bed", sep="")))
system(paste("bgzip", paste(filename_clinical, ".bed", sep=""), sep=" "))

# index ".bed"
system(paste("tabix -b 2 -e 3", paste(filename, ".bed.gz", sep=""), sep=" "))
system(paste("tabix -b 2 -e 3", paste(filename_clinical, ".bed.gz", sep=""), sep=" "))

system(paste("rm header.txt", args[["exons"]], args[["introns"]], sep=" "))

#______________________________________
# print output files
#______________________________________
cat(paste("                        ", "Output files: ", paste(filename, ".bed.gz", sep=""), " (and tabix)","\n",
"                                      ", paste(filename_clinical, ".bed.gz", sep=""), " (and tabix)", "\n",
"                                      ", "  - coverage: ", paste(filename_clinical, "_coverage.bed", sep=""), "\n", sep=""))

#______________________________________
# END LOG for each file
#______________________________________
sink()
sink()
sink()

