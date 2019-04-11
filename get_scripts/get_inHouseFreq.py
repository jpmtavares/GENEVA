#!/usr/bin/env python

import sys
import argparse
import vcf

parser = argparse.ArgumentParser(description='Get the number of samples with variantes (homozy and hetero) and the allele frequence to construct the final info')
parser.add_argument('-vcf', '--vcf_file', help='Ouptut vcf from the vcf-merge', required=True)
parser.add_argument('-outname', '--out_file', help='Name of the output file', required=True)
args = parser.parse_args()

def getAllInfo(input_vcf, out_file):

        list_var=[]
        missing_count=0
        hete_count=0
        homo_count=0
        vcf_reader = vcf.Reader(open(input_vcf, 'r'))
        col_num=len(vcf_reader.samples)
        for record in vcf_reader:
                for sample in record.samples:
                        if sample['GT']==".": missing_count=missing_count+1
                        if sample['GT']=="0/1": hete_count=hete_count+1
                        if sample['GT']=="1/1": homo_count=homo_count+1
                af=(float(2*homo_count+hete_count))/(col_num*2)
		freq="{0:.17f}".format(af)
                n=hete_count+homo_count
		chrom=record.CHROM
		pos=record.POS
		ref=record.REF
		aux=str(record.ALT)
		alt=aux.replace("[","")
		alt=alt.replace("]","")
		if "," in alt:
			i=len(alt.split(","))
			j=0
			while j<=i-1: 
				out_file.write(str(chrom) + "\t" + str(pos) + "\t" + str(ref) + "\t" + str(alt.split(",")[j]) + "\t" + str(n) + "\t" +  str(homo_count) + "\t" + str(freq) + "\n")
				j+=1
		else: out_file.write(str(chrom) + "\t" + str(pos) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(n) + "\t" +  str(homo_count) + "\t" + str(freq) + "\n")
                missing_count=0
                hete_count=0
                homo_count=0

if __name__=="__main__":

        inputFile=args.vcf_file
        outFile=args.out_file
	output=open(outFile,"w")
	getAllInfo(inputFile, output)
