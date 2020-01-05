######################################################
# This takes in a bed file without a header and outputs  
# a bed file with the longest annotations.
# TJE 2019 04 02.
#
#
#
#
######################################################

import sys
from collections import defaultdict as dd

inFile = open(sys.argv[1],'r') #input reflat, no header
outFile = open(sys.argv[2],'w+') #output refflat, no header

lineDict = dd(tuple)

#create a dictionary of distances and lines, keeping only sucessive entries
#if they have the same accession and longer distances.
for idx, line in enumerate(inFile):  
	symbol = line.split("\t")[0]
	exonStarts = line.split("\t")[9].split(",")[:-1]#the last entry is blank
	exonEnds = line.split("\t")[10].strip().split(",")[:-1] 
	strand = line.split("\t")[3]

	if strand == "+" : distance = sum([int(exonStarts[i]) - int(exonEnds[i]) \
		for i in range(0,len(exonEnds))])
	elif strand == "-": distance = sum([int(exonEnds[i]) - int(exonStarts[i]) \
		for i in range(0,len(exonEnds))])
	else: 
		print("Strand Error")
		sys.exit()

	if symbol in lineDict:
		if lineDict[symbol][1] > distance: 
			lineDict[symbol] = (line,distance)
		else: continue
	else: lineDict[symbol] = (line,distance)

#write the longest lines. 
for symbol, linehash in lineDict.items():
	outFile.write(linehash[0])

inFile.close()
outFile.close()