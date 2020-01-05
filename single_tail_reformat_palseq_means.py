from collections import defaultdict as dd
import sys
import numpy as np

file=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')
cutoff=int(sys.argv[3])

#The of the output is a file containing 252 columns, as follows:
#Accession code \t tags_of_tail_length_0 \t tags_of_tail_length_1 \t tags_of_tail_length_2...
#this is rounding Alex's single tag files. 

tag_dict=dd(list) #lambda function of list with 251 values

for line in file:
	gene=line.split("\t")[0]
	st_list=line.split("\t")[1].strip().split(",")
	if len(st_list)>=cutoff:
		for tail in st_list[:-1]:
			tag_dict[gene].append(float(tail))

for key, val in tag_dict.iteritems():
	line=str(key)+"\t"+str(np.mean(val))+"\n"
	output.write(line)
file.close()