from collections import defaultdict as dd
import sys

file=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')
cutoff=int(sys.argv[3])


#The of the output is a file containing 252 columns, as follows:
#Accession code \t tags_of_tail_length_0 \t tags_of_tail_length_1 \t tags_of_tail_length_2...

tag_dict=dd(lambda: [0]*251) #lambda function of list with 251 values

for line in file:
	tag_dict[line.split("\t")[0]][int(float(line.split("\t")[2]))]+=1

for key, val in tag_dict.iteritems():
	if sum(val)>=cutoff:
		line=str(key)+"\t"+"\t".join(map(str,val))+"\n"
		output.write(line)
file.close()
