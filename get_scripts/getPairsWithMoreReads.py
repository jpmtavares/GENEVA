#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:19:18 2019

@author: brigidameireles
"""

import argparse

parser = argparse.ArgumentParser(description='Get the number of samples with variantes (homozy and hetero) and the allele frequence to construct the final info')
parser.add_argument('-list', '--list', help='File with the number of reads of each pair', required=True)
args = parser.parse_args()

def getPairWithMoreReads(flist):
    l_reads=[]
    
    with open(flist, "r") as f:
        for line in f:
            a=line.rstrip().split("\n")
            l_reads.append(a)
    maxindex=l_reads.index(max(l_reads))
    print(maxindex)
                
if __name__=="__main__":
    
    inputf=args.list
    getPairWithMoreReads(inputf)