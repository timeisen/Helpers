#filter fastq
#this script removes reads from the fastq file that have all Bs in the QS

import sys
from Bio import SeqIO

output=open(sys.argv[2],'w+')
fastq_parser = SeqIO.parse(sys.argv[1], "fastq")
   
def my_filter(records):
    for rec in records:
    	if rec.id.split(";")[1]=="1": yield rec
    	
SeqIO.write(my_filter(fastq_parser), output, "fastq")

output.close()