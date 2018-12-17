
# GENEVA - GErmliNE Variant Analysis
A repository for annotating, interpreting, reporting and visualizing germline SNPs and small indels in exome data.

## annotations
* [RefSeq - clinical transcripts](https://github.com/jpmtavares/GENEVA/blob/master/annotations/RefSeq_clinical_transcripts.txt)

  List of RefSeq transcripts used in Human Gene Mutation Database (HGMD). This file contains information about: `HGNC_symbol`,  `ENSGene`, `ENSTranscript`, `refSeq_mRNA` and `refSeq_protein`.
 
## get_scripts
* [get_clinvar.R](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_clinvar.R)

  Get ClinVar GRCh37 vcf file from ClinVar. This script runs with default values for downloading ClinVar directory, and [vcfanno](https://github.com/brentp/vcfanno) configuration file directory. These parameters can be chosen:
  
  usage: `get_clinvar.R <ClinVar download directory> <vcfanno config file directory>`

* [get_inHousefreq.py](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_inHouseFreq.py)

  Get the number of samples in a VCF file with a given variant (either in homozigosity or heterozigosity) and the corresponding allele frequency. The output is writen in a sorted file containg: `Chr`, `Position`, `Ref`, `Alt`, `Number of samples with variant`, `Number of homozygous samples`, `Allele Frequency`.
  
  usage: `get_inHouseFreq.py -vcf <multisample VCF file> -outname <output name>`

* [get_vcfmerge2freq.sh](https://github.com/jpmtavares/GENEVA/blob/master/get_scripts/get_vcfmerge2freq.sh)

  Get multisample VCF file with [vcf-merge](http://vcftools.sourceforge.net/perl_module.html#vcf-merge), split file by chromosome, and calculates allele frequency with `get_inHousefreq.py`. This script runs with default values for every parameter:
  
    Parameter | Default value
    --- | ---
    `-v` | current directory
    `-d` | same as `-v`
    `-o` | _inHouse_freq_ _
    `-p` | 24

  
  usage: `get_vcfmerge2freq.sh [-h] [-v <path to VCF files>] [-d <path to output>] [-o <output name>] [-p <number of process to run>]`
