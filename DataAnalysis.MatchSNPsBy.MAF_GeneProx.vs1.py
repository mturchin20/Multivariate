#!/usr/bin/python

import sys
import re
import gzip
import random
from argparse import ArgumentParser

file1 = None
file2 = None

###File1 Ex: /mnt/lustre/home/mturchin20/Data/GlobalLipids2010/GlobalLipids2010.Table1.ChrBP.noNA.txt
#File1 Ex: /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2010/Vs2/bmass.GlobalLipids2010.Vs2.MergedDataSources.vs1.AnnotatedForPerms.NoPad.20kbPad.vs1.GWASSNPs.txt.gz
#File2 Ex: /mnt/lustre/home/mturchin20/Lab_Stuff/StephensLab/Multivariate/GlobalLipids2010/Vs2/bmass.GlobalLipids2010.Vs2.MergedDataSources.vs1.AnnotatedForPerms.NoPad.20kbPad.vs1.txt.gz

#Argument handling and parsing
#Parsing arguments
parser = ArgumentParser(add_help=False)

#Required arguments
required = parser.add_argument_group('required arguments:')
required.add_argument("--file1", dest="file1", help="location of file1", required=True, metavar="FILE1")
required.add_argument("--file2", dest="file2", help="location of file2", required=True, metavar="FILE2")
required.add_argument("--N", dest="N", help="number of permutations", required=True, metavar="N")
required.add_argument("--outfile1", dest="outfile1", help="location of base name outfile1", required=True, metavar="OUTFILE1")
required.add_argument("--seed", dest="seed", help="seed value for permutations", required=True, metavar="seed")

#Optional arguments
optional = parser.add_argument_group('optional arguments:')
optional.add_argument("-h", "--help", help="show this help message and exit", action="help")

args = parser.parse_args()

#print(args.file1)


#Main script

#Setting seed
random.seed(int(args.seed))

#Initializing global variables
PerPermSNPUsedRecord = []

if args.file1 == "-":
	file1 = sys.stdin
elif re.search('gz$', args.file1):
	file1 = gzip.open(args.file1, 'rb')	
else:
	file1 = open(args.file1, 'r')

for line1 in file1:
	line1 = line1.rstrip().split()	
	
	MatchedSNPs1 = []

	MAFLowerBound = int(line1[3]) - .05
	MAFUpperBound = int(line1[3]) + .05
	20kbProximityFlag = int(line1[6])
	NumberWithin250kbLowerBound = int(line1[9]) - .1*int(line1[9])
	NumberWithin250kbUpperBound = int(line1[9]) + .1*int(line1[9])

	if args.file2 == "-":
		file2 = sys.stdin
	elif re.search('gz$', args.file2):
		file2 = gzip.open(args.file2, 'rb')	
	else:
		file2 = open(args.file2, 'r')

	for line2 in file2:
		line2 = line2.rstrip().split()	

#		if int(line2[3]) >= MAFLowerBound and int(line2[3]) <= MAFUpperBound and int(line2[6]) == 20kbProximityFlag:
#		if int(line2[3]) >= MAFLowerBound and int(line2[3]) <= MAFUpperBound and int(line2[6]) == 20kbProximityFlag and int(line2[9]) >= NumberWithin250kbLowerBound and int(line2[9]) <= NumberWithin250kbLowerBound:
		if int(line2[3]) >= MAFLowerBound and int(line2[3]) <= MAFUpperBound:
			if int(line2[6]) == 20kbProximityFlag:
				if int(line2[9]) >= NumberWithin250kbLowerBound and int(line2[9]) <= NumberWithin250kbLowerBound:
					MatchedSNPs1.append(line2)

	file2.close()
	
	if len(MatchedSNPs1) <= int(args.N):
		sys.exit('Error1a -- number of available matched SNPs is less than the desired number of permutations for the original seed SNP. Error occured on: ' + ",".join(map(str, line1)))

	for i in range(args.N):
		outfile1.fullName = args.outfile1 + "Permutation" + args.N + ".txt.gz"
		with gzip.open(outfile1.fullName, 'ab') as outfile1:
			if len(PerPermSNPUsedRecord) != i+1:
				PerPermSNPUsedRecord.append({})
			PermSNPChosen = random.choice(MatchedSNPs1)
			WhileLoopCheck = 0
			WhileLoopCount = 0
			while WhileLoopCount < 100000:
				if str(PermSNPChosen[0]) not in PerPermSNPUsedRecord[i]:
					break
				else:
					PermSNPChosen = random.choice(MatchedSNPs1)
					WhileLoopCount += 1
			if WhileLoopCount == 100000:
				sys.exit('Error1b -- While loop for checking whether there are unused, matched SNPs left in the matched SNP pool reached 100000 iterations, indicating an error & possible inifinite loop occuring. Error occured on: ' + ",".join(map(str, line1)))
			outfile1.write("\t".join(map(str, PermSNPChosen)) + "\t" + str(len(MatchedSNPs1)) + "\t" + ",".join(map(str, line1)))
			PerPermSNPUsedRecord[i][str(PermSNPChosen[0])] = 1

file1.close()

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

