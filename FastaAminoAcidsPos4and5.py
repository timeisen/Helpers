#Translates the 4th and 5th amino acid of a CDS from each feature of a FASTA file. 

import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO
import re

output=open(sys.argv[2],'w+')
output.write("accession\taa4and5\n")

for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	NtAA4and5 = seq_record.seq[12:18].translate()	
	output.write(seq_record.id + "\t" + str(NtAA4and5) + "\n")

output.close()