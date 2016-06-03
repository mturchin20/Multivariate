#!/usr/bin/python

import sys
import re
import gzip
from argparse import ArgumentParser
import scipy.stats as st

file1 = None
file2 = None
hash1 = {}
hash2 = {}
hash3 = {}

#File1 Ex: /mnt/lustre/home/mturchin20/Data/ICBP2011/ICBP2011_AllPhenos_results.wHapMap21.vs3.formatted.GWASHits.MarkerChrBP.txt
#File2 Ex: /mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_AllPhenos_results.wHapMap21.vs3.formatted.txt.gz

#Argument handling and parsing
#Parsing arguments
parser = ArgumentParser(add_help=False)

#Required arguments
required = parser.add_argument_group('required arguments:')
required.add_argument("--file1", dest="file1", help="location of file1", required=True, metavar="FILE1")
required.add_argument("--file2", dest="file2", help="location of file2", required=True, metavar="FILE2")

#Optional arguments
optional = parser.add_argument_group('optional arguments:')
optional.add_argument("-h", "--help", help="show this help message and exit", action="help")

args = parser.parse_args()

#print(args.file1)

#Main script

#rs425277         1       2059032
#rs9434723        1       9214869
#rs10779751       1       11206923
#rs2284746        1       17179262
#rs12137162       1       19635983


if args.file1 == "-":
	file1 = sys.stdin
elif re.search('gz$', args.file1):
	file1 = gzip.open(args.file1, 'rb')	
else:
	file1 = open(args.file1, 'r')

for line1 in file1:
	line1 = line1.rstrip().split()	

	ChrBP1 = line1[1] + "_" + line1[2]

#	hash1[ChrBP1] = 1 
	hash1[line1[0]] = 1 

	Start1 = int(line1[2]) - 500000
	End1 = int(line1[2]) + 500000
#	Start1 = int(line1[2]) - 50
#	End1 = int(line1[2]) + 50

#	sys.stderr.write("happening1\n")

	for i in range(Start1, End1+1):
#		sys.stderr.write("hapening2: " + str(i) + "\n")
		ChrBP3 = line1[1] + "_" + str(i)
#		sys.stderr.write("hapening2: " + ChrBP3 + "\n")
#		hash3[ChrBP3] = 1
		hash3[ChrBP3] = line1

file1.close()

sys.stderr.write(args.file2 + "\n")

#[  mturchin20@spudling26  ~]$zcat /mnt/gluster/data/external_private_supp/ICBP2011/PhenoGenotypeFiles/RootStudyConsentSet_phs000585.ICBP_GWAS.v1.p1.c1.DS-CVD-IRB/AnalysisFiles/release_3_2013/data/ICBP2011_AllPhenos_results.wHapMap21.vs3.formatted.txt.gz | head -n 10
#rs4747841       10      10000135        0.476287607628189       0.785038816998804       69517.6795539446 0.0704022760997935     69513.3012699446 0.3178740130878        29131.34632     0.212811078518439 73598.0115078986
#rs4749917       10      10000265        0.47620042838397        0.784054526872054       69520.9811124215 0.0703893633932794     69516.6050754215 0.317795542112825      29130.0584025   0.212435684622182 73600.0444929245
#rs737656        10      100002729       0.363093029725855       0.539718338872868       69175.8984124676 0.939886621635876      69166.5072934676 0.461154294086731      29045.6727915   0.399477221138652 73234.9133678245
#rs737657        10      100002880       0.35788089133122        0.642748731863813       68875.554415    0.964980888514893       68866.414 0.531983750598191     28397.4596955   0.542742730514193       72942.3795295

if args.file2 == "-":
	file2 = sys.stdin
elif re.search('gz$', args.file2):
	file2 = gzip.open(args.file2, 'rb')	
else:
	file2 = open(args.file2, 'r')
	
hash4 = {}

for line2 in file2:
	line2 = line2.rstrip().split()	

	annot = "NA3"
	ChrBP2 = str(line2[1]) + "_" + str(line2[2]) 

	if line2[0] in hash4:
		sys.stderr.write("Error2a -- rsID present more than once (line2: " + ",".join(map(str, line2)) + ")\n")
	else:
		hash4[line2[0]] = 1
		
	if line2[0] in hash1:
		annot = 1
		line2.append(annot)
		line2.extend(hash3[ChrBP2])
	elif ChrBP2 in hash3:
		annot = 2
		line2.append(annot)
		line2.extend(hash3[ChrBP2])
	else:
		annot = 0

	line2.append(annot)
	print "\t".join(map(str, line2))


file2.close()

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

