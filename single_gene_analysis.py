import numpy as np
import sys
from collections import defaultdict as dd
#this script takes in a singletails output file and makes a similar one with just the gene name and tail length for genes that 
#have more than "cutoff" tags

file=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')
cutoff=100

long_list=dd(list)
for line in file:
	long_list[line.split("\t")[0]].append(float(line.split("\t")[3].strip()))

for key, val in long_list.iteritems():
	if len(val)>=cutoff:
		for tail_length in val:
			line=str(key)+"\t"+str(tail_length)+"\n"
			output.write(line)
file.close()
