#filter fastq
#this script returns records that match an input file id

import sys
from Bio import SeqIO

filter_file=open(sys.argv[1],'r')
fastq_parser = SeqIO.parse(sys.argv[2], "fastq")
output=open(sys.argv[3],'w+')

filter_set=set()
for line in filter_file:
	filter_set.add(line.strip())

def checkEqual2(iterator):
   return len(set(iterator)) <= 1

def my_filter(records):
    for rec in records:
    	id=":".join(rec.id.split(":")[2:5]).split("#")[0]
    	if id in filter_set: yield rec

SeqIO.write(my_filter(fastq_parser), output, "fastq")

output.close()