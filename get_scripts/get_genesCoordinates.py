#!/usr/bin/python

import argparse
import gzip

parser = argparse.ArgumentParser(
    description='Get the genes position, first and last position')
parser.add_argument('-genes_file', '--refseq', help='RefSeq file, sorted in bed format', required=True)
parser.add_argument('-outname', '--out_file', help='Name of the output file', required=True)
args = parser.parse_args()


def createFileWithGenes(infile, outname):
    dic_genes = {}
    outfile=open(outname, "w")

    with gzip.open(infile, "rb") as bed:
        next(bed)
        for line in bed:
            aux = line.rstrip().split("\t")
            id = aux[5]
            if "_" not in aux[0]:
                if id in dic_genes.keys():
                    info = dic_genes.get(id, " ")
                    dic_genes[id] = [info[0], info[1], aux[2], info[3]]
                else:
                    dic_genes[aux[5]] = [aux[0], aux[1], aux[2], aux[6]]


    outfile.write("#Chro"+"\t"+"Start"+"\t"+"End"+"\t"+"Gene"+"\t"+"Gene_alt"+"\n")
    for key, values in dic_genes.items():
        outfile.write(values[0]+"\t"+values[1]+"\t"+values[2]+"\t"+key+"\t"+values[3]+"\n")

    outfile.close()
    bed.close()

if __name__ == "__main__":
    f_genes = args.refseq
    outname = args.out_file

createFileWithGenes(f_genes, outname)
