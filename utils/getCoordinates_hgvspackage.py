# -*- coding: utf-8 -*-
"""
Created on Fri Mar  1 11:58:04 2019

@author: brigidameireles
"""
import argparse
import hgvs.variantmapper
import hgvs.parser
import hgvs.assemblymapper
import hgvs.dataproviders.uta
import re

parser = argparse.ArgumentParser(description='Use HGVS package to find the Chromosome and the Position')
parser.add_argument('-in', '--input', help='Input example: \'NM_001637.3:c.1582G>A\'', required=True)
args = parser.parse_args()

if __name__ == "__main__":

    hgvs_c = args.input
    d_chromo={'NC_000001':1, 'NC_000002': 2, 'NC_000003':3, 'NC_000004':4, 'NC_000005':5, 'NC_000006':6, 'NC_000007':7, 'NC_000008': 8, 'NC_000009':9, 'NC_000010':10, 'NC_000011':11, 'NC_000012': 12, 'NC_000013':13, 'NC_000014':14, 'NC_000015':15, 'NC_000016':16, 'NC_000017':17, 'NC_000018':18, 'NC_000019':19, 'NC_000020':20,'NC_000021':21,'NC_000022':22,'NC_000023':'X', 'NC_000024':'Y'} 
    h=hgvs.parser.Parser()
    hdp = hgvs.dataproviders.uta.connect()
    variantmapper = hgvs.variantmapper.VariantMapper(hdp)
    variantmapper = hgvs.assemblymapper.AssemblyMapper(hdp, assembly_name="GRCh37")
    
    var_c=h.parse_hgvs_variant(hgvs_c)
    var_p=variantmapper.c_to_g(var_c)
    refseq_chromosome=str(var_p).split('.')[0]
    chromo=d_chromo[refseq_chromosome]
    position=re.findall('\d+', (str(var_p).split(':')[1]))
    pos_a = (str(var_p).split('>')[1])
    pos_b=(str(var_p).split('>')[0][-1])
    print("Chomosome: " + str(chromo))
    print("Position: " + position[0])
    print("Ref: " + pos_b + "    Alt: "+pos_a)