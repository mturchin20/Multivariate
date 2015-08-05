#!/usr/bin/python

import sys
import re
import gzip
from argparse import ArgumentParser

file1 = None
file2 = None
hash1 = {}

#File1 Ex: /mnt/lustre/home/mturchin20/Data/dbSNP/All_20150415.vcf.gz
#File2 Ex: /mnt/gluster/data/external_public_supp/GIANT2010/GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt.gz

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

#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
#rs10 a c 0.0333 0.708 80566
#rs1000000 g a 0.6333 0.506 123865
#rs10000010 c t 0.425 0.736 123827
#rs10000012 c g 0.8083 0.042 123809

if args.file1 == "-":
	file1 = sys.stdin
elif re.search('gz$', args.file1):
	file1 = gzip.open(args.file1, 'rb')	
else:
	file1 = open(args.file1, 'r')

for line1 in file1:
	line1 = line1.rstrip().split()	

	if line1[0] in hash1:
		sys.stderr.write("Error1a -- there are multiple entries of a rsID (" + str(line1[0]) + ")\n")
	else:
		hash1[line1[0]] = line1

file1.close()

#1       10250   rs199706086     A       C       .       .       RS=199706086;RSPOS=10250;dbSNPBuildID=137;SSR=0;SAO=0;VP=0x050000020005000002000100;WGT=1;VC=SNV;R5;ASP
#1       10254   rs140194106     TA      T       .       .       RS=140194106;RSPOS=10255;dbSNPBuildID=134;SSR=0;SAO=0;VP=0x050000020005000002000200;WGT=1;VC=DIV;R5;ASP
#1       10257   rs111200574     A       C       .       .       RS=111200574;RSPOS=10257;dbSNPBuildID=132;SSR=0;SAO=0;VP=0x050100020005000102000100;WGT=1;VC

if args.file2 == "-":
	file2 = sys.stdin
elif re.search('gz$', args.file2):
	file2 = gzip.open(args.file2, 'rb')	
else:
	file2 = open(args.file2, 'r')

for line2 in file2:
	line2 = line2.rstrip().split()	

	if not re.search('^#', line2[0]):
	
		ChrBP2 = str(line2[0]) + "_" + str(line2[1])

		if line2[2] in hash1:
#			print "\t".join(hash1[line2[2]]) + "\t" + ChrBP2
			hash1[line2[2]].append(ChrBP2)

file2.close()

for entry1 in hash1:
	if len(hash1[entry1]) == 6:
		hash1[entry1].append("NA")
	
	print "\t".join(hash1[entry1])

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

