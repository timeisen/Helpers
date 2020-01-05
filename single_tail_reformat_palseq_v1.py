from collections import defaultdict as dd
import sys
import numpy as np

file=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')
cutoff=int(sys.argv[3])

#The of the output is a file containing 252 columns, as follows:
#Accession code \t tags_of_tail_length_0 \t tags_of_tail_length_1 \t tags_of_tail_length_2...
#this is rounding Alex's single tag files. 

tag_dict=dd(lambda: [0]*251) #lambda function of list with 251 values

for line in file:
	gene=line.split("\t")[0]
	st_list=line.split("\t")[1].strip().split(",")
	if len(st_list)>=cutoff:
		for tail in st_list[:-1]:
			if int(round(float(tail)))>=250: tag_dict[gene][250]+=1
			elif int(round(float(tail)))<=0: tag_dict[gene][0]+=1
			else: tag_dict[gene][int(round(float(tail)))]+=1

for key, val in tag_dict.iteritems():
	if sum(val)>=cutoff:
		line=str(key)+"\t"+"\t".join(map(str,val))+"\n"
		output.write(line)
file.close()