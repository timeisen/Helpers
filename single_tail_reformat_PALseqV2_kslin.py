from collections import defaultdict as dd
import sys
#PYTHON 3!!
file=open(sys.argv[1],'r')
output=open(sys.argv[2],'w+')
cutoff=int(sys.argv[3])

#The of the output is a file containing 252 columns, as follows:
#Accession code \t tags_of_tail_length_0 \t tags_of_tail_length_1 \t tags_of_tail_length_2...
next(file)
header = "\t".join(["tl"+str(i) for i in range(0,251)]) #make a fancy header row.
output.write("accession\t"+header+"\n")
tag_dict=dd(lambda: [0]*251) #lambda function of list with 251 values

for line in file:
	accession = line.split("\t")[6].strip()
	tag_dict[accession][int(line.split("\t")[2])]+=1

for key, val in tag_dict.items():
	if sum(val)>=cutoff:
		line=str(key)+"\t"+"\t".join(map(str,val))+"\n"
		output.write(line)
file.close()
