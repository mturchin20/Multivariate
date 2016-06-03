#!/usr/bin/python

import sys
import re
import gzip
from argparse import ArgumentParser

file1 = None
file2 = None
hash1 = {}

##File1 Ex: /mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release22/Frequencies/allele_freqs_chrAll_CEU_r22_nr.b36.txt.gz
#File2 Ex: /mnt/gluster/data/external_public_supp/HaemgenRBC2012/

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

##HapMap allele_frequency file
#rs# chrom pos strand build center protLSID assayLSID panelLSID QC_code refallele refallele_freq refallele_count otherallele otherallele_freq otherallele_count totalcount
#rs11511647 chr10 62765 + ncbi_b36 sanger urn:lsid:illumina.hapmap.org:Protocol:Golden_Gate_1.0.0:1 urn:lsid:sanger.hapmap.org:Assay:4310385:1 urn:lsid:dcc.hapmap.org:Panel:CEPH-30-trios:1 QC+ T 0.178 21 C 0.822 97 118
#rs4880608 chr10 83299 + ncbi_b36 affymetrix urn:LSID:bcm.hapmap.org:Protocol:genotype_0002:1 urn:LSID:bcm.hapmap.org:Assay:319692:1 urn:lsid:dcc.hapmap.org:Panel:CEPH-30-trios:1 QC+ G 1 116 A 0 0 116
#rs12218882 chr10 84172 + ncbi_b36 perlegen urn:lsid:perlegen.hapmap.org:Protocol:Genotyping_1.0.0:2 urn:lsid:perlegen.hapmap.org:Assay:25770.5651698:1 urn:lsid:dcc.hapmap.org:Panel:CEPH-30-trios:1 QC+ G 0.967 116 A 0.033 4 120
#rs10904045 chr10 84426 + ncbi_b36 perlegen urn:lsid:perlegen.hapmap.org:Protocol:Genotyping_1.0.0:2 urn:lsid:perlegen.hapmap.org:Assay:25770.7352214:1 urn:lsid:dcc.hapmap.org:Panel:CEPH-30-trios:1 QC+ C 0.593 70 T 0.407 48 118

Chrs = ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"]

for chr in Chrs:

#	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release21/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r21_nr.txt.gz"
#	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release22/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r22_nr.b36.txt.gz"
#	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release23/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r23_nr.b36_fwd.txt.gz"
#	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release24/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r24_nr.b36_fwd.txt.gz"
#	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release27/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r27_nr.b36_fwd.txt.gz"
	args.file1 = "/mnt/lustre/home/mturchin20/Data/HapMap/Releases/Release28/Frequencies/allele_freqs_chr" + str(chr) + "_CEU_r28_nr.b36_fwd.txt.gz"

	if args.file1 == "-":
		file1 = sys.stdin
	elif re.search('gz$', args.file1):
		file1 = gzip.open(args.file1, 'rb')	
	else:
		file1 = open(args.file1, 'r')

	for line1 in file1:
		line1 = line1.rstrip().split()	

		if not re.search('^#', line1[0]):
	#		ChrBP1 = str(line1[0]) + "_" + str(line1[1])
	#
	#		if line1[2] in hash1:
	#			sys.stderr.write("Error1a -- there are multiple entries of a rsID (" + str(line1[2]) + ")\n")
	#		else:
	#			hash1[line1[2]] = ChrBP1

	#		ChrBP1 = line1[1].split('hr')[1].split('_')[0] + "_" + str((int(line1[2]) + int(round(float((int(line1[3]) - int(line1[2])))/2))))
	
			AFInfo1 = str(line1[1]) + "_" + str(line1[2]) + "_" + str(line1[10]) + "_" + str(line1[11]) + "_" + str(line1[13]) + "_" + str(line1[14])

			if line1[0] in hash1:
				hash1[line1[0]].append(AFInfo1)			
	#			sys.stderr.write("Error1a -- there are multiple entries of a rsID (" + str(line1[0]) + ")\n")
			else:
				hash1[line1[0]] = [AFInfo1]

	file1.close()

#MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
#rs10 a c 0.0333 0.708 80566
#rs1000000 g a 0.6333 0.506 123865
#rs10000010 c t 0.425 0.736 123827
#rs10000012 c g 0.8083 0.042 123809


if args.file2 == "-":
	file2 = sys.stdin
elif re.search('gz$', args.file2):
	file2 = gzip.open(args.file2, 'rb')	
else:
	file2 = open(args.file2, 'r')

for line2 in file2:
	line2 = line2.rstrip().split()	

	if re.search('^rs\d+', line2[0]):
		if line2[0] in hash1:
			line2.append(",".join(map(str, hash1[line2[0]])))
		else:
			line2.append("NA")
		
		print "\t".join(line2)
	else:
		line2.append("ChrBPAFInfo")
		print "\t".join(line2)

file2.close()


def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

