# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 16:56:38 2019

@author: brigidameireles
"""

#!/usr/bin/env python
import argparse
import datetime
import re
import sys

parser = argparse.ArgumentParser(description='Get the sample name (input: forward and reverse raw reads). If the sample name is ausent of the raw fastq write the name of the sample in the option -n. If the tool used in preprocessed is changed is necessary change the output of this script')
parser.add_argument('-f', '--pair1', help='Forward fastq', required=True)
parser.add_argument('-r', '--pair2', help='Reverse fastq', required=True)
parser.add_argument('-n', '--name', help='Sample name', required=False)
args = parser.parse_args()

def checkformat(file_format):
    form=""    
    if "fastq" in file_format: form="fastq"
    if "fq" in file_format: form="fq"
    else:
        if "fq" not in file_format or "fastq" not in file_format: 
            print("ERROR! " + "The input files are in a format that is not fastq or fq or are not compressed. Check the input, correct and do a new run.")
            sys.exit()        
    return form
    
def checkpairs (f_raw, r_raw):
    outname = "log_"
    outf=open(outname, "w")
    outf.write("[ "+str(datetime.datetime.now().replace(microsecond=0))+" ] "+"STEP1: Checking the input."+"\n")

    try: form=checkformat(f_raw.split(".")[-2])
    except IndexError:         
        print("ERROR! "+"Not a valid input, check the fastq files.")
        sys.exit()
    outf.write("                        " + "The "+form+" are compressed and in correct format."+"\n")
    whichpair_f=f_raw.split(".")[-3][-6:]
    whichpair_r=r_raw.split(".")[-3][-6:]
    if re.match("^R1.*", whichpair_f):
        if re.match("^R2.*", whichpair_r): outf.write("                        " + "The input is pair 1 and pair 2, respectively."+"\n")
        else: 
            print("ERROR! "+"Check the second pair, it is not concordant with the first pair.") 
            sys.exit()        
    if re.match(".*[._]1.*", whichpair_f) or re.match("^1", whichpair_f):
        if re.match(".*[._]2.*", whichpair_r) or re.match("^2", whichpair_r): outf.write("                        " + "The input is pair 1 and pair 2, respectively."+"\n")
        else: 
            print("ERROR! "+"Check the second pair, it is not concordant with the first pair.") 
            sys.exit()
    if not re.match("^R1.*", whichpair_f) and not re.match(".*[._]1.*", whichpair_f) and not re.match("^1", whichpair_f):
            print("ERROR! "+"Check the first pair, it is not correct.") 
            sys.exit()
    return outf
        
def getName_SampleWithoutName (sample_name, outf):
    outf.write("[ "+str(datetime.datetime.now().replace(microsecond=0))+" ] "+"Preparing the reads for the preprocessing."+"\n")
    firstpair=(sample_name.upper()+"_1.fastp.fq.gz")
    secondpair=(sample_name.upper()+"_2.fastp.fq.gz")
    return firstpair, secondpair, sample_name
    
def getName_check (pair1, outf):
    outf.write("[ "+str(datetime.datetime.now().replace(microsecond=0))+" ] "+"Preparing the reads for the preprocessing."+"\n")
    firstpair=""
    secondpair=""
    if pair1.startswith("Sample"): #Sample68006_HM_FDHE19H000077-1a-A60-A67_HWYWFCCXY_L3_1.fq.gz
        name=pair1.rstrip().split("_")[0:2]
        sample_name=name[0].replace("Sample","")+"-"+name[1]
        firstpair=(sample_name+"_1.fastp.fq.gz")
        secondpair=(sample_name+"_2.fastp.fq.gz")
    if len(str(pair1.rstrip().split("_")[0]))==6 and len(pair1.rstrip().split("_")[1])==2: #h66887_MH_DHT02688-A4-A61_HT3TNCCXY_L1_2.fq.gz
        if re.match("^[a-zA-Z]+.*", pair1): #starts with any letter
            name=pair1.rstrip().split("_")[0:2]
            first=str(name[0])[0]
            sample_name=name[0].replace(first,"")+"-"+name[1]
            firstpair=(sample_name+"_1.fastp.fq.gz")
            secondpair=(sample_name+"_2.fastp.fq.gz")
    if len(str(pair1.rstrip().split("_")[0]))==2 and len(pair1.rstrip().split("_")[1])==5: #PS_67859_DHE17173-1-A96-A56_H7Y5CDSXX_L3_1.fq.gz
        if str(pair1.rstrip().split("_")[0]).isalpha() and pair1.rstrip().split("_")[1].isdigit():
            sample_name=pair1.rstrip().split("_")[1]+"-"+pair1.rstrip().split("_")[0]
            firstpair=(sample_name+"_1.fastp.fq.gz")
            secondpair=(sample_name+"_2.fastp.fq.gz")
    return firstpair, secondpair, sample_name
        
if __name__ == "__main__":
    firstpair_raw= args.pair1
    secondpair_raw = args.pair2    
    
    outf=checkpairs(firstpair_raw, secondpair_raw)
    if args.name:
        sample_name=args.name
        forward,reverse, sample_name=getName_SampleWithoutName(sample_name, outf)
    else: forward, reverse, sample_name=getName_check(firstpair_raw, outf)
    outf.write("                        ==>First pair:   " + forward +"\n")
    outf.write("                        ==>Second pair:  " + reverse+"\n")
    outf.write("[ "+str(datetime.datetime.now().replace(microsecond=0))+" ] "+"STEP2: Starting the preprocessing with FASTP tool."+"\n")
    print(sample_name+","+forward+","+reverse)
    
