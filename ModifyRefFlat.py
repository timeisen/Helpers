################################################################################
# ModifyRefFlat.py
# Timothy J Eisen
# 2020 01 05
# Takes in a refFlat file and a list of 3' end annotations from PAL-seq data
# Modifies the txEnd coordinates of the refFlat file, keeping genes that are not
# in the PAL-seq data unchanged. 
# Only uses the most-used 3p end.
# Format of input annotation: accession\tchrom\tstrand\t3pEndPos\ttagFraction\n
################################################################################


import sys, pdb
from collections import defaultdict as dd

RefFlat  = open(sys.argv[1],'r')
PalAnno = open(sys.argv[2],'r')
RefFlatOut   = open(sys.argv[3],'w+')

#Skip the header
next(PalAnno)

#Iterate through the pal-seq annotations, creating a dictionary of all 3p ends
# for each gene
AnnoDict = dd(lambda: ([],[]))
for line in PalAnno:
	#Parse the line
	accession       = line.split("\t")[0]
	chrom		    = line.split("\t")[1]
	strand          = line.split("\t")[2]
	txEnd           = int(line.split("\t")[3])
	tagFrac         = float(line.split("\t")[4].strip())

	AnnoDict[accession][0].append(txEnd)
	AnnoDict[accession][1].append(tagFrac)

#Iterate through the pal-seq annotation dictionary, keeping only max tag 3' end
TxEndDict = {}
for accession, tagTuple in AnnoDict.items():
	maxTag = max(tagTuple[1])
	maxPos = [i for i, j in enumerate(tagTuple[1]) if j == maxTag][0]
	TxEndDict[accession] = tagTuple[0][maxPos]

#iterate through RefFlat
for line in RefFlat:

	pos = True

	#Parse the line.
	strand        = line.split("\t")[3]
	chrom         = line.split("\t")[2]
	accession     = line.split("\t")[1]
	NumberOfExons = int(line.split("\t")[8])

	txstart       = int(line.split("\t")[4])
	txEnd         = int(line.split("\t")[5])
	cdsStart      = int(line.split("\t")[6])
	cdsEnd        = int(line.split("\t")[7])
	exonStarts    = line.split("\t")[9].split(",")[:-1] #remove last element
	exonEnds      = line.strip().split("\t")[10].split(",")[:-1]

	if accession not in TxEndDict:
		RefFlatOut.write(line)
		continue

	if strand == ('-'): pos = False

	#Change tx end to the correct coordinate above.
	#Change the last exon to end at this coordinate
	#Remove downstream exons
	#If it is before CdsEnd, throw a warning with the accession code and make CdsEnd = TxEnd.

#Close the files
RefFlat.close()
PalAnno.close()
RefFlatOut.close()
