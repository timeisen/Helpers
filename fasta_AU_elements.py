#Calculate the AU element percentage from each feature of a FASTA file. 
#Definitition and code from Floor and Doudna, 2016.  

import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO
import re

output=open(sys.argv[2],'w+')
output.write("accession\tAU_PentamerCount\tnumAUElements" + \
			"\taurichFraction\tlongestAUElement\n")

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	AUcounts = seq_record.seq.count("G") + seq_record.seq.count("C")


	# calculate: 
	# 1) count AUUUA
	# 2) number of AU-stretches
	# 3) percentage of UTR covered by AREs #A COVERAGE DEFINITION
	# 4) longest A/U stretch 
	auPentamerCount = seq_record.seq.count("ATTTA") 

	aurich_re = re.compile(r'[AT]{5,}')  # RE to find more than 5 a/u in a row 
	au_elements = aurich_re.findall(str(seq_record.seq))

	if (au_elements):
		numAUElements = len(au_elements)
		aurichFraction = len("".join(au_elements))/float(len(seq_record.seq))
		longestAUElement = max(map(len,au_elements))

	else: 
		numAUElements = aurichFraction = longestAUElement = 0

	output.write(seq_record.id + "\t" + \
		str(auPentamerCount) + "\t" + \
		str(numAUElements) + "\t" + \
		str(aurichFraction) + "\t" +\
		str(longestAUElement) + "\n")

output.close()