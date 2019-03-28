#!/usr/bin/Rscript

#______________________________________________
# libraries
#______________________________________________
library(stringr, quietly = T)
library(plyr, quietly = T)
library(dplyr, quietly = T, warn.conflicts = F)
library(magrittr, quietly = T)
library(ape, quietly = T)
library(valr, quietly = T)
library(biomaRt, quietly = T)
library(lubridate)
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
#______________________________________
# START LOG for each file
#______________________________________
sink(file = "log_grch37.refseq_ensembl_lrg_hugo.txt", append = FALSE, type = c("output"),
     split = TRUE)

cat("########################################################\n")
cat("grch37.refseq_ensembl_lrg_hugo.R\n")
cat("########################################################\n")
cat(paste("[", now("GMT"), "]", "Creating GRCh37 RefSeq transcripts annotation file, with Ensembl, LRG, HUGO genes and clinical transcript information", "\n"))

#______________________________________
# print input files
#______________________________________
cat(paste("                       ", "Input files:", "GRCh37_latest_genomic.gff.gz", "(from: https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml)","\n", 
"                                   ", "Ensembl.grch37vs38_genes.txt", "(from: GENEVA)", "\n",
"                                   ", "Ensembl BIOMART", "(R package)", "\n",
"                                   ", "LRG", "(from: ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_GRCh37.bed)", "\n",
"                                   ", "grch37.clin.manual.refseq_ensembl.txt", "(from: GENEVA)", "\n",
"                                   ", "HUGO genes", "(from: ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt)", "\n", sep=" "))

################################################################
# STEP1: Getting RefSeq GRCh37 annotation from NCBI
################################################################
#https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml

system("wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz")

refseq_gff<-read.gff("GRCh37_latest_genomic.gff.gz",
                     na.strings = c(".", "?"), GFF3 = TRUE) %>%
  #filter(type=="mRNA" | type=="CDS") %>%
  cbind(data.frame(str_split_fixed(.$attributes, ";", 10)))

# 1) get gene names (they can be at columns 6, 7, 8 and 9)
gene_name6<-unlist(str_extract(refseq_gff$X6, "(?<=gene\\=)(.*?)$"))
gene_name7<-unlist(str_extract(refseq_gff$X7, "(?<=gene\\=)(.*?)$"))
gene_name8<-unlist(str_extract(refseq_gff$X8, "(?<=gene\\=)(.*?)$"))
gene_name9<-unlist(str_extract(refseq_gff$X9, "(?<=gene\\=)(.*?)$"))

gene_name<-ifelse(is.na(gene_name8), gene_name6, gene_name8) %>%
  ifelse(is.na(.), gene_name7, .) %>%
  ifelse(is.na(.), gene_name9, .)

rm(gene_name6, gene_name7, gene_name8, gene_name9)

# 2) get NM_ , NP_  and NR_ transcripts (they are at column 4)
refseq_nm_np_nr<-unlist(str_extract(refseq_gff$X4, "(?<=Name\\=)(.*?)$"))

# 3) join two tables
refseq_df<-data.frame(gene_name, refseq_nm_np_nr) %>%
  unique %>%
  filter(!is.na(refseq_nm_np_nr)) %>% #remove genes without NM_ or NP_
  filter(!grepl("YP", refseq_nm_np_nr)) #remove genes with YP_

RefSeq<-data.frame(HGNC_symbol=refseq_df$gene_name[grep("NM_", refseq_df$refseq_nm_np_nr)],
                   refSeq_mRNA=refseq_df$refseq_nm_np_nr[grep("NM_", refseq_df$refseq_nm_np_nr)],
                   refSeq_protein=refseq_df$refseq_nm_np_nr[grep("NP_", refseq_df$refseq_nm_np_nr)],
                   refSeq_mRNA_noVersion=str_extract(refseq_df$refseq_nm_np_nr[grep("NM_", refseq_df$refseq_nm_np_nr)],
                                                     "(.*?)(?=\\.)"),
                   refSeq_protein_noVersion=str_extract(refseq_df$refseq_nm_np_nr[grep("NP_", refseq_df$refseq_nm_np_nr)],
                                                        "(.*?)(?=\\.)"))

RefSeq_nr<-data.frame(HGNC_symbol=refseq_df$gene_name[grep("NR_", refseq_df$refseq_nm_np_nr)],
                      refSeq_mRNA=refseq_df$refseq_nm_np_nr[grep("NR_", refseq_df$refseq_nm_np_nr)],
                      refSeq_protein=NA,
                      refSeq_mRNA_noVersion=str_extract(refseq_df$refseq_nm_np_nr[grep("NR_", refseq_df$refseq_nm_np_nr)],
                                                        "(.*?)(?=\\.)"),
                      refSeq_protein_noVersion=NA)

RefSeq<-rbind(RefSeq, RefSeq_nr)

rm(refseq_df, gene_name, refseq_nm_np_nr, RefSeq_nr)
system("rm GRCh37_latest_genomic.gff.gz")

################################################################
# STEP2: Getting GRCh37 vs GRCh38 list of genes that changed their names
################################################################
system("wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/annotations/Ensembl.grch37vs38_genes.txt")

grch37vs38_table<-read.delim("Ensembl.grch37vs38_genes.txt", header=T)

final38<-inner_join(RefSeq, grch37vs38_table, by=c(HGNC_symbol="GRCh38"))
names(final38)[ncol(final38)]<-"HGNC_alternative_symbol"
final37<-inner_join(RefSeq, grch37vs38_table, by=c(HGNC_symbol="GRCh37"))
names(final37)[ncol(final37)]<-"HGNC_alternative_symbol"

grch37vs38<-rbind(final38, final37) %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol) %>%
  unique %>%
  left_join(RefSeq, .) %>%
  mutate(HGNC_alternative_symbol = ifelse(is.na(HGNC_alternative_symbol),
                                          HGNC_symbol,
                                          as.character(HGNC_alternative_symbol)))
rm(final37, final38, grch37vs38_table)
system("rm Ensembl.grch37vs38_genes.txt")
################################################################
# STEP3: Getting Ensembl transcript information from BioMart
################################################################
#GRCh37
#Get gene names with NR and NM separately
nr_genes<-unique(c(grch37vs38$HGNC_symbol[grep("NR_", grch37vs38$refSeq_mRNA)],
                   grch37vs38$HGNC_alternative_symbol[grep("NR_", grch37vs38$refSeq_mRNA)]))

nm_genes<-unique(c(grch37vs38$HGNC_symbol[grep("NM_", grch37vs38$refSeq_mRNA)],
                   grch37vs38$HGNC_alternative_symbol[grep("NM_", grch37vs38$refSeq_mRNA)]))

#Get Ensembl BioMart
ensembl<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                 path="/biomart/martservice" , dataset="hsapiens_gene_ensembl")
########
#  NM  #
########
Ensembl_info_nm<-getBM(attributes = c('chromosome_name','hgnc_symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'refseq_mrna', 'refseq_peptide'), 
                    filters = 'hgnc_symbol', 
                    values = nm_genes, 
                    mart = ensembl) %>%
  filter(chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                                "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")) %>%
  dplyr::select(-chromosome_name)
names(Ensembl_info_nm)<-c("HGNC_symbol_ensembl", "ENSGene", "ENSTranscript", "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")
########
#  NP  #
########
Ensembl_info_nr<-getBM(attributes = c('chromosome_name','hgnc_symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'refseq_mrna'), 
                       filters = 'hgnc_symbol', 
                       values = nr_genes, 
                       mart = ensembl) %>%
  filter(chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                                "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")) %>%
  dplyr::select(-chromosome_name) %>%
  mutate(refSeq_protein_noVersion=NA)
names(Ensembl_info_nr)<-c("HGNC_symbol_ensembl", "ENSGene", "ENSTranscript", "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")
#####################
#  join two tables  #
#####################
Ensembl_info<-rbind(Ensembl_info_nm, Ensembl_info_nr)
# inner join HGNC_symbol_ensembl == HGNC_symbol
Ensembl_symbol<-inner_join(grch37vs38, Ensembl_info, by=c("HGNC_symbol"="HGNC_symbol_ensembl",
                                                          "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")) 
# inner join HGNC_symbol_ensembl == HGNC_alternative_symbol
Ensembl_alternative_symbol<-inner_join(grch37vs38, Ensembl_info, by=c("HGNC_alternative_symbol"="HGNC_symbol_ensembl",
                                                                      "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")) 
# add two tables
Ensembl_join<-rbind(Ensembl_symbol, Ensembl_alternative_symbol) %>%
  unique

# left join with grch37vsgrch38 in order to get entries that didn't match before
Ensembl<-left_join(grch37vs38, Ensembl_join, by=c("HGNC_symbol", "HGNC_alternative_symbol",
                                                  "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")) %>%
  unique %>%
  dplyr::select(-refSeq_mRNA.y, -refSeq_protein.y)

names(Ensembl)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
                  "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript")

rm(Ensembl_alternative_symbol, Ensembl_info, Ensembl_info_nr, Ensembl_info_nm, Ensembl_join, Ensembl_symbol)
################################################################
# STEP4: LRG IDs
################################################################
lrg_bed<-read_bed("ftp://ftp.ebi.ac.uk/pub/databases/lrgex/LRG_GRCh37.bed",
                  n_fields = 12)

# 1) get LRG IDs
lrg_id<-unlist(str_extract(lrg_bed$name, "(.*?)(?=\\()"))
# 2) get gene names, and NM_ transcripts
gene_name<-unlist(str_extract(lrg_bed$name, "(?<=\\()(.*?)(?=\\)|\\|)"))

# 3) join two tables
lrg_df<-data.frame(lrg_id, gene_name) %>%
  filter(grepl("NM_", gene_name)) # get only NM_ transcripts
names(lrg_df)<-c("LRG_id", "refSeq_mRNA")

lrg_df$refSeq_mRNA_noVersion<-str_extract(lrg_df$refSeq_mRNA, "(.*?)(?=\\.)")

# left join with Ensembl and lrg_df
LRG<-left_join(Ensembl, lrg_df, by=c("refSeq_mRNA_noVersion")) %>%
  unique %>%
  dplyr::select(-refSeq_mRNA.y)

names(LRG)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
              "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript", "LRG_id")

rm(lrg_bed, lrg_df, lrg_id)

################################################################
# STEP5: Getting clinical transcripts
################################################################
#system("wget https://raw.githubusercontent.com/jpmtavares/GENEVA/master/annotations/grch37.clin.manual.refseq_ensembl.txt")

clinical_table<-read.delim(paste(MENDEL,"clinical_manual/grch37.clin.manual.refseq_ensembl.txt",sep=""), header=T)

#left join with LRG and clinical transcripts 
# get clinical transcripts without ENSEMBL IDs
clinical_withoutEnsembl<-LRG %>%
  filter(is.na(ENSTranscript)) %>%
  left_join(., clinical_table, by=c("refSeq_mRNA_noVersion", "refSeq_protein_noVersion")) %>%
  unique %>%
  #if a gene is already in the clinical transcripts list
  #it was manually curated -> "yes_manual"
  mutate(HGNC_symbol.y=ifelse(is.na(HGNC_symbol.y), "no", "yes_manual")) %>%
  dplyr::select(-ENSGene.y, -ENSTranscript.y)

names(clinical_withoutEnsembl)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
                                  "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript",
                                  "LRG_id", "clinical_transcript")

# get clinical transcripts with ENSEMBL IDs
clinical_withEnsembl<-LRG %>%
  filter(!is.na(ENSTranscript)) %>%
  left_join(., clinical_table, by=c("refSeq_mRNA_noVersion", "refSeq_protein_noVersion",
                                    "ENSGene", "ENSTranscript")) %>%
  unique %>%
  #if a gene is already in the clinical transcripts list
  #it was manually curated -> "yes_manual"
  mutate(HGNC_symbol.y=ifelse(is.na(HGNC_symbol.y), "no", "yes_manual"))

names(clinical_withEnsembl)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
                               "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript",
                               "LRG_id", "clinical_transcript")

# join clinical transcripts with and without ENSEMBL IDs
clinical<-rbind(clinical_withEnsembl, clinical_withoutEnsembl) %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol,
                refSeq_mRNA, refSeq_protein,
                refSeq_mRNA_noVersion, refSeq_protein_noVersion,
                ENSGene, ENSTranscript,
                LRG_id, clinical_transcript) %>%
  unique

rm(clinical_table, clinical_withEnsembl, clinical_withoutEnsembl)
system("rm grch37.clin.manual.refseq_ensembl.txt")

################################################################
# STEP6: Getting HUGO genes transcripts
################################################################
system("wget ftp://ftp.ebi.ac.uk/pub/databases/genenames/new/tsv/locus_types/gene_with_protein_product.txt")

hugo_table<-read.delim("gene_with_protein_product.txt", header=T) %>%
  dplyr::select(symbol, ensembl_gene_id, refseq_accession)
names(hugo_table)<-c("HGNC_symbol", "ENSGene", "refSeq_mRNA_noVersion")

HUGO<-clinical %>%
  left_join(., hugo_table, by=c("refSeq_mRNA_noVersion", "ENSGene")) %>%
  group_by(HGNC_symbol.x) %>%
  #if a gene doesn't have clinical transcripts (all(clinical_transcript=="no"))
  #and it's described in HUGO gene -> "yes_HUGO"
  mutate(clinical_transcript=ifelse((all(clinical_transcript=="no") & !(is.na(HGNC_symbol.y))), "yes_HUGO", clinical_transcript)) %>%
  mutate(n_transcript=n()) %>%
  #if a gene doesn't have clinical transcripts
  #and has only 1 transcript, that transcript is its clinical transcript -> "yes_automatic"
  mutate(clinical_transcript=ifelse((n_transcript==1 & clinical_transcript=="no"), "yes_automatic",
                                    #if a gene has more than one transcript
                                    #and has only 1 coding transcript, that transcript is its clinical transcript -> "yes_coding"
                                    ifelse((n_transcript>1 & length(grep("NM_", refSeq_mRNA))==1), "yes_coding", clinical_transcript))) %>%
  ungroup() %>%
  #if a gene has more than one transcript
  #and has only 1 coding transcript -> "yes_coding", all the remaining NR_ transcripts are not clinical -> "yes_coding"
  mutate(clinical_transcript=ifelse((clinical_transcript=="yes_coding" & grepl("NR_", refSeq_mRNA)), "no", clinical_transcript)) %>%
  dplyr::select(-HGNC_symbol.y,-n_transcript)

names(HUGO)<-c("HGNC_symbol", "HGNC_alternative_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
               "refSeq_protein_noVersion", "ENSGene", "ENSTranscript",
               "LRG_id", "clinical_transcript")

rm(hugo_table)
system("rm gene_with_protein_product.txt")
################################################################
# STEP7: Getting warnings
################################################################
#______________________________________________
# WARNINGS...
#______________________________________________
geneswithoutNM_Ensembl<-Ensembl %>%
  group_by(HGNC_symbol) %>%
  mutate(geneswithoutNM_Ensembl=ifelse(any(is.na(ENSGene) == FALSE),
                                       "all_good", "check_this_gene")) %>%
  ungroup %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol, geneswithoutNM_Ensembl) %>%
  unique

geneswithoutclinical<-HUGO %>%
  group_by(HGNC_symbol) %>%
  mutate(geneswithoutclinical=ifelse(any(clinical_transcript %in% "no" == FALSE),
                                     "all_good", "check_this_gene")) %>%
  ungroup %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol, geneswithoutclinical) %>%
  unique

warning<-cbind(geneswithoutclinical, geneswithoutNM_Ensembl$geneswithoutNM_Ensembl)
names(warning)[ncol(warning)]<-"geneswithoutNM_Ensembl"

################################################################
# STEP8: Output files
################################################################
#______________________________________________
# write.tables
#______________________________________________
write.table(HUGO, "grch37.refseq_ensembl_lrg_hugo.txt", col.names=T,
            row.names=F, quote = F, sep="\t")
write.table(warning, "grch37.refseq_ensembl_lrg_hugo.warnings.txt", col.names=T,
            row.names=F, quote = F, sep="\t")

#______________________________________
# print output files
#______________________________________
cat(paste("                        ", "Output files: ", "grch37.refseq_ensembl_lrg_hugo.txt", "\n",
"                                      ", "grch37.refseq_ensembl_lrg_hugo.warnings.txt", "\n", sep=""))

sink()
