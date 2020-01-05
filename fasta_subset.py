#This script takes the following inputs

import sys, argparse
from collections import defaultdict as dd
from Bio import SeqIO

##Parse arguments
parser=argparse.ArgumentParser(description='subsets a fasta file based on a list of genes')
parser.add_argument('--gene_file','-g', action='store',dest="gene_file", help='A list of accession codes, one per line')
parser.add_argument('--fasta','-f', action='store',dest="fasta_file", help='Fasta file with identifiers matching the gene list')
parser.add_argument('--output','-o', action='store',dest="output_file", help='Output fasta file')
args = parser.parse_args()

gene_list=open(args.gene_file,'r')
#utr_list="/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sequences/Mouse/mm10_tpUTR_fixed.txt"
output=open(args.output_file,'w+')

k=[]

for line in gene_list:
	k.append(line.strip())
for seq_record in SeqIO.parse(args.fasta_file, "fasta"):
	if seq_record.id in k: 
		output.write(">"+str(seq_record.id)+"\n"+str(seq_record.seq)+"\n")

gene_list.close()
output.close()