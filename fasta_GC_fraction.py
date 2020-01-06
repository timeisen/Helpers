#Calculate the GC percentage from each feature of a FASTA file. 

import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO

output=open(sys.argv[2],'w+')
output.write("accession\tGC_counts\tGC_fraction\n")
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	GCcounts = seq_record.seq.count("G") + seq_record.seq.count("C")
	output.write(seq_record.id + "\t" + str(GCcounts) + \
		"\t" + str(float(GCcounts) / len(seq_record.seq)) +\
		"\n")


output.close()