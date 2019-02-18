#______________________________________________
# libraries
#______________________________________________
library(stringr, quietly = T)
library(plyr, quietly = T)
library(dplyr, quietly = T)
library(magrittr, quietly = T)
library(ape, quietly = T)
library(valr, quietly = T)
library(fuzzyjoin, quietly = T)
library(biomaRt, quietly = T)
#______________________________________________
# create RefSeq annotation NM_ and NP_ 
#______________________________________________
print("Getting RefSeq GRCh37 annotation from NCBI:")
print("ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz")
#https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml
system("wget ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh37_latest/refseq_identifiers/GRCh37_latest_genomic.gff.gz")
print("NCBI RefSeq GRCh37 retrieved.")

system("Reading GFF file...")
refseq_gff<-read.gff("GRCh37_latest_genomic.gff.gz",
                     na.strings = c(".", "?"), GFF3 = TRUE) %>%
  filter(type=="mRNA" | type=="CDS") %>%
  cbind(data.frame(str_split_fixed(.$attributes, ";", 10)))

# 1) get gene names (they can be at columns 6, 7, 8 and 9)
system("Getting Gene names and RefSeq transcripts.")
gene_name6<-unlist(str_extract(refseq_gff$X6, "(?<=gene\\=)(.*?)$"))
gene_name7<-unlist(str_extract(refseq_gff$X7, "(?<=gene\\=)(.*?)$"))
gene_name8<-unlist(str_extract(refseq_gff$X8, "(?<=gene\\=)(.*?)$"))
gene_name9<-unlist(str_extract(refseq_gff$X9, "(?<=gene\\=)(.*?)$"))

gene_name<-ifelse(is.na(gene_name8), gene_name6, gene_name8) %>%
  ifelse(is.na(.), gene_name7, .) %>%
  ifelse(is.na(.), gene_name9, .)

rm(gene_name6, gene_name7, gene_name8, gene_name9)

# 2) get NM_ and NP_ transcripts (they are at column 4)
refseq_nm_np<-unlist(str_extract(refseq_gff$X4, "(?<=Name\\=)(.*?)$"))

# 3) join two tables
refseq_df<-data.frame(gene_name, refseq_nm_np) %>%
  unique %>%
  filter(!is.na(refseq_nm_np)) %>% #remove genes without NM_ or NP_
  filter(!grepl("YP", refseq_nm_np)) #remove genes with YP_

system("Output NCBI RefSeq transcripts.")
RefSeq<-data.frame(HGNC_symbol=refseq_df$gene_name[grep("NM_", refseq_df$refseq_nm_np)],
                   refSeq_mRNA=refseq_df$refseq_nm_np[grep("NM_", refseq_df$refseq_nm_np)],
                   refSeq_protein=refseq_df$refseq_nm_np[grep("NP_", refseq_df$refseq_nm_np)],
                   refSeq_mRNA_noVersion=str_extract(refseq_df$refseq_nm_np[grep("NM_", refseq_df$refseq_nm_np)],
                                                     "(.*?)(?=\\.)"),
                   refSeq_protein_noVersion=str_extract(refseq_df$refseq_nm_np[grep("NP_", refseq_df$refseq_nm_np)],
                                                        "(.*?)(?=\\.)"))

rm(refseq_df, gene_name, refseq_nm_np)
system("rm GRCh37_latest_genomic.gff.gz")
#______________________________________________
# alternative names
#______________________________________________
system("Reading GRCh37 vs GRCh38 genes that changed their names...")
grch37vs38_table<-read.delim("grch37vs38.txt", header=T)

final38<-inner_join(RefSeq, grch37vs38_table, by=c(HGNC_symbol="GRCh38"))
names(final38)[ncol(final38)]<-"HGNC_alternative_symbol"
final37<-inner_join(RefSeq, grch37vs38_table, by=c(HGNC_symbol="GRCh37"))
names(final37)[ncol(final37)]<-"HGNC_alternative_symbol"

system("Join GRCh37vsGRCh38 list of genes with NCBI RefSeq transcripts output.")
grch37vs38<-rbind(final38, final37) %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol) %>%
  unique %>%
  left_join(RefSeq, .) %>%
  mutate(HGNC_alternative_symbol = ifelse(is.na(HGNC_alternative_symbol),
                                          HGNC_symbol,
                                          as.character(HGNC_alternative_symbol)))
rm(final37, final38)
system("rm grch37vs38.txt")
#______________________________________________
# get Ensembl IDs
#______________________________________________
print("Getting Ensembl transcript information from BioMart...")
#GRCh37
ensembl<-useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", 
                 path="/biomart/martservice" , dataset="hsapiens_gene_ensembl")
Ensembl_info<-getBM(attributes = c('chromosome_name','hgnc_symbol', 'ensembl_gene_id', 'ensembl_transcript_id', 'refseq_mrna', 'refseq_peptide'), 
                    filters = 'hgnc_symbol', 
                    values = unique(c(grch37vs38$HGNC_symbol, grch37vs38$HGNC_alternative_symbol)), 
                    mart = ensembl) %>%
  filter(chromosome_name %in% c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", 
                                "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")) %>%
  dplyr::select(-chromosome_name)
names(Ensembl_info)<-c("HGNC_symbol_ensembl", "ENSGene", "ENSTranscript", "refSeq_mRNA_noVersion", "refSeq_protein_noVersion")
print("BioMart retrieved.")

system("Join Ensembl information with NCBI RefSeq transcripts output.")
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

#______________________________________________
# get LRG IDs
#______________________________________________
system("Getting LRG transcripts...")
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
system("Join LRG transcripts with NCBI RefSeq transcripts output.")
LRG<-left_join(Ensembl, lrg_df, by=c("refSeq_mRNA_noVersion")) %>%
  unique %>%
  dplyr::select(-refSeq_mRNA.y)

names(LRG)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
              "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript", "LRG_id")
#______________________________________________
# get clinical transcripts
#______________________________________________
system("Getting clinical transcripts...")
clinical_table<-read.delim("RefSeq_clinical_transcripts.txt", header=T)

#left join with LRG and clinical transcripts 
system("Join clinical transcripts with NCBI RefSeq transcripts output.")
clinical<-left_join(LRG, clinical_table, by=c("refSeq_mRNA_noVersion", "refSeq_protein_noVersion",
                                              "ENSGene", "ENSTranscript")) %>%
  unique %>%
  mutate(HGNC_symbol.y=ifelse(is.na(HGNC_symbol.y), "no", "yes"))

names(clinical)<-c("HGNC_symbol", "refSeq_mRNA", "refSeq_protein", "refSeq_mRNA_noVersion",
                   "refSeq_protein_noVersion", "HGNC_alternative_symbol", "ENSGene", "ENSTranscript",
                   "LRG_id", "clinical_transcript")

clinical<-clinical %>% 
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol,
                refSeq_mRNA, refSeq_protein,
                refSeq_mRNA_noVersion, refSeq_protein_noVersion,
                ENSGene, ENSTranscript,
                LRG_id, clinical_transcript)

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

geneswithoutclinical<-clinical %>%
  group_by(HGNC_symbol) %>%
  mutate(geneswithoutclinical=ifelse(any(clinical_transcript %in% "no" == FALSE),
                                  "all_good", "check_this_gene")) %>%
  ungroup %>%
  dplyr::select(HGNC_symbol, HGNC_alternative_symbol, geneswithoutclinical) %>%
  unique

warning<-cbind(geneswithoutclinical, geneswithoutNM_Ensembl$geneswithoutNM_Ensembl)
names(warning)[ncol(warning)]<-"geneswithoutNM_Ensembl"

#______________________________________________
# write.tables
#______________________________________________
write.table(clinical, "RefSeqGRCh37_Ensembl_LRG_clinical.txt", col.names=T, row.names=F, quote = F, sep="\t")
write.table(warning, "RefSeqGRCh37_Ensembl_LRG_clinical.checklist.txt", col.names=T, row.names=F, quote = F, sep="\t")
