#sum up tags for pal-seq 3p UTR comparisons

from collections import defaultdict as dd
import sys

bedfile=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')

tag_dict=dd(lambda: 0) #lambda function of list with 251 values

for line in bedfile:

	strand = line.split("\t")[5]
	chr = line.split("\t")[0]
	left_coordinate = line.split("\t")[1]
	right_coordinate = line.split("\t")[2]

	if int(right_coordinate)-int(left_coordinate)!=50: continue #include this line to require that reads are 50nt exactly. 

	if strand == "+":
		hash = right_coordinate+strand+chr
 	if strand == "-":
		hash = left_coordinate+strand+chr

	tag_dict[hash]+=1

for key, val in tag_dict.iteritems():
	line=key+"\t"+str(val)+"\n"
	output.write(line)
bedfile.close()
output.close()