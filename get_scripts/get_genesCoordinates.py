import argparse

parser = argparse.ArgumentParser(
    description='Cross the information from HGMD and from HGVS, completing the HGMD with HGVS nomenclature')
parser.add_argument('-genes_file', '--refseq', help='RefSeq file, sorted in bed format', required=True)
parser.add_argument('-outname', '--out_file', help='Name of the output file', required=True)
args = parser.parse_args()


def createFileWithGenes(infile, outname):
    dic_genes = {}
    outfile=open(outname, "w")

    with open(infile) as bed:
        for line in bed:
            aux = line.rstrip().split("\t")
            id = aux[5]
            if id in dic_genes.keys():
                info = dic_genes.get(id, " ")
                dic_genes[id] = [info[0], info[1], aux[2]]
            else:
                dic_genes[aux[5]] = [aux[0], aux[1], aux[2]]


    outfile.write("#Chro"+"\t"+"Start"+"\t"+"End"+"\t"+"Gene"+"\n")
    for key, values in dic_genes.items():
        outfile.write(values[0]+"\t"+values[1]+"\t"+values[2]+"\t"+key+"\n")

    outfile.close()
    bed.close()

if __name__ == "__main__":
    f_genes = args.refseq
    outname = args.out_file

    createFileWithGenes(f_genes, outname)
