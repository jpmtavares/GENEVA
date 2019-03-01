import sys
import argparse
import re
from seqseek import Chromosome, BUILD37

parser = argparse.ArgumentParser(description='Get the HGVS nomenclature for each variant.')
parser.add_argument('-vcf', '--in_vcf', help='Ouptut from VEP, version 95', required=True)
parser.add_argument('-outname', '--out_file', help='Name of the output file', required=True)
args = parser.parse_args()

def getAlleleRef(chromo, posi):

        if chromo=="X" or chromo=="Y":
                if chromo=="X": #give an error, the function didn't recognize the 'X'
                        allele=Chromosome('X').sequence((int(posi)-1),int(posi))
                else:
                        allele=Chromosome('Y').sequence((int(posi)-1),int(posi))
        else:
                allele=Chromosome(chromo).sequence((int(posi)-1),int(posi))
		return allele

def getOutput(vcfin, vcfout):
	lvars=[]
	hgvsp="-"
	hgvsc="-"
	np="-"
	nm="-"
	warn=False #if need, create a warning file
	
	outerro=vcfout+".warning.vcf"
	outfile=vcfout+".vcf"
	firstLine=True
	output=open(outfile,"w")
	vcf=open(vcfin, "r")
	for line in vcf:
		if not line.startswith("#") and "HGVS" in line:
			aux=line.rstrip().split("\t")[13]
			rs_id=line.rstrip().split("\t")[0]
			chro_pos=line.rstrip().split("\t")[1]
			chro=chro_pos.split(":")[0]
			pos=chro_pos.split(":")[1]
			if "-" in pos: pos=pos.split("-")[0]
			alt=line.rstrip().split("\t")[2]
			inf=aux.split(";")
			gref=inf[2].split("=")[1]
			uref=inf[3].split("=")[1]
			if gref != uref:
				if firstLine:
					outwarn=open(outerro, "w")
					firstLine=False
					outwarn.write(line)
				warn=True
			else:
				if "-" in gref: gref=getAlleleRef(chro, pos)
				for i in range(len(inf)) :
					if 'HGVSc' in inf[i]:
						if ":" in inf[i]:
							nm=inf[i].split(":")[0].replace("HGVSc=","")
							hgvsc=inf[i].split(":")[1]
						else:	nm=inf[i].replace("HGVSc=","")
					if 'HGVSp' in inf[i]:
						if ":" in inf[i]:
							np=inf[i].split(":")[0].replace("HGVSp=","")
							hgvsp=inf[i].split(":")[1]
						else:	np=inf[i].replace("HGVSp=","")
				if "%" in hgvsp: hgvsp=hgvsp.replace("%3D","=")
				output.write(chro+"\t"+pos+"\t"+gref+"\t"+alt+"\t"+nm+"\t"+hgvsc+"\t"+np+"\t"+hgvsp+"\t"+rs_id+"\n")
	output.close()
	if warn:
		outwarn.close()
				
if __name__=="__main__":
    inputFile=args.in_vcf
    outname=args.out_file

    getOutput(inputFile, outname)
