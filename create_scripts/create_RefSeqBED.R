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
      create_RefSeqBED.R [--help] [--exons==<path to UCSC BED file with exons>] [--introns==<path to UCSC BED file with introns>] [--RefSeq==<path to file with custom RefSeq annotation>]
      
      -- program to create BED file with RefSeq clinical transcripts along with ranking for exons and introns --
      
      where:
      --exons                      path to UCSC BED file with exons
      --introns                    path to UCSC BED file with introns
      --RefSeq                     path to file with custom RefSeq annotation\n\n")
  
  q(save="no")
}

# SETUP
#________________________________________________________
print("Reading input files...")
exons<-read.delim(args[["exons"]], header= F)
introns<-read.delim(args[["introns"]], header= F)
RefSeq<-read.delim(args[["RefSeq"]], header = TRUE)
filename<-"RefSeqGRCh37"
filename_clinical<-"RefSeqGRCh37_clinical"
print("Done.")
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
print("Correct Region order in minus Strand...")
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
print("Add clinical transcript infomation...")
bed_clinical<-bed %>%
  filter(clinical_transcript=="yes") %>%
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
dir.create("RefSeq_annotation/", showWarnings = TRUE, recursive = FALSE)
print("Output files in RefSeq_annotation/...")

write.table(bed, paste("RefSeq_annotation/", filename, ".bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t")
write.table(bed_clinical, paste("RefSeq_annotation/", filename_clinical, ".bed", sep=""), col.names=F, row.names=F, quote=F, sep="\t")

# add header to final files
write.table("#Chr\tStart\tEnd\tRegion\tStrand\tHGNC_symbol\tHGNC_alternative_symbol\tENSGene\tENSTranscript\trefSeq_mRNA\trefSeq_protein\trefSeq_mRNA_noVersion\trefSeq_protein_noVersion\tLRG_id",
            "RefSeq_annotation/header.txt", col.names=F, row.names=F, quote=F, sep="\t")
system(paste("cat RefSeq_annotation/header.txt", paste("RefSeq_annotation/", filename, ".bed", sep=""), ">",
             paste("RefSeq_annotation/", filename, "_hdr.bed", sep=""), sep=" "))
system(paste("cat RefSeq_annotation/header.txt", paste("RefSeq_annotation/", filename_clinical, ".bed", sep=""), ">",
             paste("RefSeq_annotation/", filename_clinical, "_hdr.bed", sep=""), sep=" "))
system("rm RefSeq_annotation/header.txt")
################################################################
# 6) sort BED file, create BED for coverage analysis and indexed it
################################################################
#print("Only the files with clinical transcripts are being sorted and tabixed...")
#print("Uncomment script to get full file RefSeqGRCh37.bed sorted and tabixed.")
#system(paste("less", paste("RefSeq_annotation/", filename, "_hdr.bed", sep="") , "|", "body sort -k1,1 -k2,2n", "-", ">", paste("RefSeq_annotation/", filename, "_hdr_sort.bed", sep=""), sep=" "))
#system(paste("less", paste("RefSeq_annotation/", filename_clinical, "_hdr.bed", sep="") , "|", "body sort -k1,1 -k2,2n", "-", ">", paste("RefSeq_annotation/", filename_clinical, "_hdr_sort.bed", sep=""), sep=" "))

# BED for coverage analysis: Chr    Start    End    Gene,Region,Strand,ENSGene,ENSTranscript,RefSeqNM
#system(paste("awk 'NF{NF-=1};1'", paste("RefSeq_annotation/", filename, "_hdr_sort.bed", sep=""), "| uniq | awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"$6\",\"$4\",\"$5\",\"$8\",\"$9\",\"$10}' | sed -r 's/\\s+/\\t/g' >", paste("RefSeq_annotation/", filename, "_hdr_coverage.bed", sep=""), sep=" "))
#system(paste("awk 'NF{NF-=1};1'", paste("RefSeq_annotation/", filename_clinical, "_hdr_sort.bed", sep=""), "| uniq | awk '{ print $1\"\t\"$2\"\t\"$3\"\t\"$6\",\"$4\",\"$5\",\"$8\",\"$9\",\"$10}' | sed -r 's/\\s+/\\t/g' >", paste("RefSeq_annotation/", filename_clinical, "_hdr_coverage.bed", sep=""), sep=" "))

# index "_sort.bed"
#system(paste("bgzip", paste("RefSeq_annotation/", filename, "_sort.bed", sep=""), sep=" "))
#system(paste("bgzip", paste("RefSeq_annotation/", filename_clinical, "_hdr_sort.bed", sep=""), sep=" "))

#system(paste("tabix -b 2 -e 3", paste("RefSeq_annotation/", filename, "_sort.bed.gz", sep=""), sep=" "))
#system(paste("tabix -b 2 -e 3", paste("RefSeq_annotation/", filename_clinical, "_hdr_sort.bed.gz", sep=""), sep=" "))
print("Finished! Bye.")
