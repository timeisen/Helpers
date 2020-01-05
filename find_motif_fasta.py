import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO

##Parse arguments
parser=argparse.ArgumentParser(description='Writes identifiers of fasta entries that contain a motif')
parser.add_argument('--fasta','-f', action='store',dest="fasta_file", help='Fasta file with identifiers matching the gene list')
parser.add_argument('--output','-o', action='store',dest="output_file", help='Output gene list')
args = parser.parse_args()

#utr_list="/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sequences/Mouse/mm10_tpUTR_fixed.txt"
output=open(args.output_file,'w+')

k=[]
motif_list=['TTTTTTTT','TTTTTGTT','TTTTTCTT']

for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
	for motif in motif_list:
		if motif in seq_record.seq: output.write(str(seq_record.id)+"\n")


output.close()