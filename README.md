
# GENEVA - GErmliNE Variant Analysis
A repository for annotating, interpreting, reporting and visualizing germline SNPs and small indels in exome data.

## annotations
* [RefSeq](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt)

  List of RefSeq GRCh37 transcripts integrated with Ensembl, LRG and [clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt). This file was created with [create_RefSeqGRCh37_Ensembl_LRG_clinical.R](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/create_RefSeqGRCh37_Ensembl_LRG_clinical.R) and contains information about: `HGNC_symbol`,  `HGNC_alternative_symbol`, `refSeq_mRNA`, `refSeq_protein`, `refSeq_mRNA_noVersion`, `refSeq_protein_noVersion`, `ENSGene`, `ENSTranscript`, `LRG_id`, `clinical_transcript`.
  
* [RefSeq - clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt)

  List of RefSeq transcripts used in Human Gene Mutation Database (HGMD). This file was manually curated and contains information about: `HGNC_symbol`,  `ENSGene`, `ENSTranscript`, `refSeq_mRNA` and `refSeq_protein`.
  
* [GRCh37vs38](https://github.com/jpmtavares/GENEVA/blob/master/annotations/Ensembl.grch37vs38_genes.txt)
  
  List of genes that changed their names between genome versions GRCh37 and GRCh38. This list was retrieved from Ensembl. 
  
## create_scripts
* [create_RefSeqBED.R](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/create_RefSeqBED.R)

  Creates NCBI RefSeq BED file with clinical transcripts information. This script accepts RefSeq BED files (hg19) from UCSC Table Browser `--exons=="path/to/exons.bed"` and `--introns=="path/to/introns.bed"`, as well as a BED file with [RefSeq, Ensembl, LRG and clinical information](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt) `--RefSeq==path/to/RefSeqGRCh37_Ensembl_LRG_clinical.txt`. The output file contains information about: `Chr`, `Start 1-based`, `End`, `Rank.Exons.Introns`, `Strand`, `HGNC_symbol`, `HGNC_alternative_symbol`,  `ENSGene`, `ENSTranscript`, `refSeq_mRNA` and `refSeq_protein`, `refSeq_mRNA_noVersion`, `refSeq_protein_noVersion`, `LRG_id`.
  
  **usage**: `create_RefSeqBED.R --exons=="path/to/exons.bed" --introns=="path/to/introns.bed" --RefSeq=="path/to/RefSeqGRCh37_Ensembl_LRG_clinical.txtt"`
 
 * [create_RefSeqGRCh37_Ensembl_LRG_clinical.R](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/create_RefSeqGRCh37_Ensembl_LRG_clinical.R)
 
   Creates [RefSeqGRCh37_Ensembl_LRG_clinical.txt](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt). This script downloads NCBI RefSeq GRCh37 from [here](https://www.ncbi.nlm.nih.gov/projects/genome/guide/human/index.shtml), integrates [GRCh37vs38](https://github.com/jpmtavares/GENEVA/blob/master/annotations/grch37vs38.txt), gets Ensembl Gene and Transcripts IDs from BioMart, gets LRG GRCh37 transcripts from [here](https://www.lrg-sequence.org/data/), and integrates manually curated [clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt) from HGMD. It outputs all these information in `RefSeqGRCh37_Ensembl_LRG_clinical.txt` and writes a list of genes without any clinical transcript defined and/or a list of genes without RefSeq_mRNA in `RefSeqGRCh37_Ensembl_LRG_clinical.checklist.txt`. 
   
    **usage**: `create_RefSeqGRCh37_Ensembl_LRG_clinical.R`
 
## get_scripts
* [get_clinvar.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_clinvar.R)

  Gets ClinVar GRCh37 vcf file from ClinVar. This script runs with default values for downloading ClinVar directory, and [vcfanno](https://github.com/brentp/vcfanno) configuration file directory. These parameters can be chosen:
  
  **usage**: `get_clinvar.R <ClinVar download directory> <vcfanno config file directory>`
  
 * [get_genesCoordinates.py](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_genesCoordinates.py) 

    Gets the first and last position of a gene. The output is written in a sorted filed by coordinates: `Chr`, `First Position`, `Last Position`, `Gene`. 
    
    **usage**: `get_genesCoordinates.py [-h] -genes_file REFSEQ -outname OUT_FILE`

* [get_hgvsnomenclature.py](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_hgvsnomenclature.py)  

    Gets the variants that have [HGVS](http://www.hgvs.org/) information (cDNA and protein) from VEP (v.95) output (vcf) with the human genome (hg19). The output is a vcf with: `Chr`, `First Position`, `Last Position`, `RefSeq_mRNA`, `HGVS c.`, `RefSeq_prot`, `HGVS p.`. 
    
    **usage**: `get_hgvsnomenclature.py [-h] -vcf VCFFile -outname OUT_FILE`

* [get_inHousefreq.py](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_inHouseFreq.py)

  Gets the number of samples in a VCF file with a given variant (either in homozigosity or heterozigosity) and the corresponding allele frequency. The output is writen in a sorted file containg: `Chr`, `Position`, `Ref`, `Alt`, `Number of samples with variant`, `Number of homozygous samples`, `Allele Frequency`.
  
  **usage**: `get_inHouseFreq.py -vcf <multisample VCF file> -outname <output name>`
  
* [get_mutationtaster.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_mutationtaster.R)

  Gets [MutationTaster](http://www.mutationtaster.org/) prediction scores. The output doesn't get transcript information into consideration, presenting unique values for a certain variant. This script accepts an input file with `Chr`, `Position`, `Ref`, `Alt`.
  
  **usage**: `get_mutationtaster.R --variants=="path/to/variants_file.txt"`
              can be parallelized with: `find <path/to/variants_files*> | xargs -n1 -P 20 -I {} sh -c 'echo {} && ./get_mutationtaster.R --variants=={}'`

* [get_vcfmerge2freq.sh](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_vcfmerge2freq.sh)

  Gets multi-sample VCF file with [vcf-merge](http://vcftools.sourceforge.net/perl_module.html#vcf-merge), split file by chromosome, and calculates allele frequency with `get_inHousefreq.py`. This script runs with default values for every parameter:
  
    Parameter | Default value
    --- | ---
    `-v` | current directory
    `-d` | same as `-v`
    `-o` | _inHouse_freq_ _
    `-p` | 24

  
  **usage**: `get_vcfmerge2freq.sh [-h] [-v <path to VCF files>] [-d <path to output>] [-o <output name>] [-p <number of process to run>]`
  
* [get_UMDpredictor.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_UMDpredictor.R)

  Gets [UMD-predictor](http://umd-predictor.eu/) scores for a list of Ensembl transcripts IDs. The output is written in a sorted file containg: `Chr`, `Position`, `Ref`, `Alt`, `HGVS_c`, `HGVS_p`, `HGNC_symbol`, `ENSTranscript`, `UMD_pred`, `UMD_score`

  **usage**: `get_UMDpredictor.R [--help] --ENSTranscripts=="path/to/ENSTranscripts.txt"`
  
  ## pipelines
* [Allele Frequency](https://github.com/jpmtavares/GENEVA/blob/master/pipelines/GENEVA_AlleleFrequency.sh)

  Downloads and runs `get_vcfmerge2freq.sh` and `get_inHouseFreq.py` in present work directory, which will be used as default directory for getting VCF samples files and output Allele Frequency file.
  
  **usage**: `GENEVA_AlleleFrequency.sh`
  
* [RefSeq BED annotation files](https://github.com/jpmtavares/GENEVA/blob/master/pipelines/GENEVA_RefSeqBED.sh)

  Downloads RefSeq BED files (hg19) from UCSC Table Browser `exons.bed` and `introns.bed`, as well as, [RefSeqGRCh37_Ensembl_LRG_clinical.txt](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeqGRCh37_Ensembl_LRG_clinical.txt) created with [create_RefSeqGRCh37_Ensembl_LRG_clinical.R](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/create_RefSeqGRCh37_Ensembl_LRG_clinical.R). It also downloads and runs `create_RefSeqBED.R` and `get_genesCoordinates.py` in present work directory. This pipeline outputs the following files in `RefSeq_annotation/` directory:

1) **clinical/** (uncomment script to write the same outputs for complete BED)
  
  
  * `RefSeqGRCh37_clinical_coordinates.txt`, file with start and end of each gene in `RefSeqGRCh37_clinical_sort.bed.gz` useful for tabix in posterior analyses.
  * `RefSeqGRCh37_clinical_coverage.bed`, input file to use in [coverage analysis](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/createCoverageDoc.py).
  * `RefSeqGRCh37_clinical_hdr_sort.bed.gz` and `RefSeqGRCh37_clinical_hdr_sort.bed.gz.tbi`, sorted and indexed BED file with clinical transcripts
  
  **usage**: `GENEVA_RefSeqBED.sh`
  
  ## utils
* [checkDiffNameGene.py](https://github.com/jpmtavares/GENEVA/blob/master/utils/checkDiffNameGene.py) 

  This script crosses the name of the genes in both versions of human genome (GRCh37 and GRCh38) and write a file with the genes that has changed.
  
  **usage**:`checkGenesNamesgrch37vsgrch38.py [-h] -refSeq REF_FILE -twoVersions BOTH_FILE -outname OUT_FILE`

* [vcf2table.R](https://github.com/jpmtavares/GENEVA/blob/master/utils/vcf2table.R)

  Converts a VCF file into a TSV with the corresponding header. The output file is written in the same directory as the input.
  
  **usage**: `vcf2table.R --vcf=="path/to/input_file.vcf" [optional: --header=="CHROM, POS, ID, REF, ALT"]`
