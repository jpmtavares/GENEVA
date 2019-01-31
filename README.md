
# GENEVA - GErmliNE Variant Analysis
A repository for annotating, interpreting, reporting and visualizing germline SNPs and small indels in exome data.

## annotations
* [RefSeq - clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt)

  List of RefSeq transcripts used in Human Gene Mutation Database (HGMD). This file contains information about: `HGNC_symbol`,  `ENSGene`, `ENSTranscript`, `refSeq_mRNA` and `refSeq_protein`.
  
## create_scripts
* [create_RefSeq.R](https://github.com/jpmtavares/GENEVA/blob/master/create_scripts/create_RefSeq.R)

  Create NCBI RefSeq BED file with clinical transcripts information. This script accepts RefSeq BED files (hg19) from UCSC Table Browser `--exons=="path/to/exons.bed"` and `--introns=="path/to/introns.bed"`, as well as a BED file with [clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt) `--clinical==path/to/clinical_transcripts.txt`. The output file contains information about: `Chr`, `Start 1-based`, `End`, `Rank.Exons.Introns`, `Strand`, `HGNC_symbol`,  `ENSGene`, `ENSTranscript`, `refSeq_mRNA` and `refSeq_protein`.
  
  **usage**: `create_RefSeq.R --exons=="path/to/exons.bed" --introns=="path/to/introns.bed" --clinical=="path/to/clinical_transcripts.txt"`
 
## get_scripts
* [get_clinvar.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_clinvar.R)

  Get ClinVar GRCh37 vcf file from ClinVar. This script runs with default values for downloading ClinVar directory, and [vcfanno](https://github.com/brentp/vcfanno) configuration file directory. These parameters can be chosen:
  
  **usage**: `get_clinvar.R <ClinVar download directory> <vcfanno config file directory>`

* [get_inHousefreq.py](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_inHouseFreq.py)

  Get the number of samples in a VCF file with a given variant (either in homozigosity or heterozigosity) and the corresponding allele frequency. The output is writen in a sorted file containg: `Chr`, `Position`, `Ref`, `Alt`, `Number of samples with variant`, `Number of homozygous samples`, `Allele Frequency`.
  
  **usage**: `get_inHouseFreq.py -vcf <multisample VCF file> -outname <output name>`

* [get_vcfmerge2freq.sh](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_vcfmerge2freq.sh)

  Get multisample VCF file with [vcf-merge](http://vcftools.sourceforge.net/perl_module.html#vcf-merge), split file by chromosome, and calculates allele frequency with `get_inHousefreq.py`. This script runs with default values for every parameter:
  
    Parameter | Default value
    --- | ---
    `-v` | current directory
    `-d` | same as `-v`
    `-o` | _inHouse_freq_ _
    `-p` | 24

  
  **usage**: `get_vcfmerge2freq.sh [-h] [-v <path to VCF files>] [-d <path to output>] [-o <output name>] [-p <number of process to run>]`
  
* [get_UMDpredictor.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_UMDpredictor.R)

  Get [UMD-predictor](http://umd-predictor.eu/) scores for a list of Ensembl transcripts IDs. The output is written in a sorted file containg: `Chr`, `Position`, `Ref`, `Alt`, `HGVS_c`, `HGVS_p`, `HGNC_symbol`, `ENSTranscript`, `UMD_pred`, `UMD_score`

  **usage**: `get_UMDpredictor.R [--help] --ENSTranscripts=="path/to/ENSTranscripts.txt"`
  
  ## pipelines
* [Allele Frequency](https://github.com/jpmtavares/GENEVA/blob/master/pipelines/GENEVA_AlleleFrequency.sh)

  Download and run `get_vcfmerge2freq.sh` and `get_inHouseFreq.py` in present work directory, which will be used as default directory for getting VCF samples files and output Allele Frequency file.
  
  **usage**: `GENEVA_AlleleFrequency.sh`
