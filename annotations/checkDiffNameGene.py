import argparse

parser = argparse.ArgumentParser(description='Check the names of the genes in the two versions and output if the name change (grch37 and grch38).')
parser.add_argument('-refSeq', '--ref_file', help='File with the transcripts of reference', required=True)
parser.add_argument('-twoVersions', '--both_file', help='File with the genes info from both versions of the genome', required=True)
parser.add_argument('-outname', '--out_file', help='Name of the output file', required=True)
args = parser.parse_args()

def getNameByVersion(genesNameF):
    dic_cross={}

    with open(genesNameF) as verFile:
        next(verFile)
        for line in verFile:
            aux=line.rstrip().split("\t")
            if aux[0] not in dic_cross:
                dic_cross[aux[0]]=(aux[1], aux[2])
    verFile.close()
    return dic_cross

def getRefName(refSeqF):
    dic_ref={}

    with open(refSeqF) as inputFile:
        for line in inputFile:
            aux=line.rstrip().split("\t")
            if aux[1] not in dic_ref:
                dic_ref[aux[1]]=aux[0]
    inputFile.close()
    return dic_ref

def checkInNameChange(dic_NamesRef, dic_NamesBothVers, outname):

    outf = open(outname, "w")
    outf.write('#ENSGene'+"\t"+'GRCh38'+"\t"+'GRCh37'+"\n")
    for key in dic_NamesRef.keys():
        if key in dic_NamesBothVers.keys():
            outf.write(str(key)+"\t"+dic_NamesBothVers[key][0]+"\t"+dic_NamesBothVers[key][1]+"\n")

if __name__ == "__main__":
    refseqFile = args.ref_file
    bothFile = args.both_file
    outname = args.out_file

    dic_NamesBothVers=getNameByVersion(bothFile)
    dic_NamesRef=getRefName(refseqFile)
    checkInNameChange(dic_NamesRef, dic_NamesBothVers, outname)
