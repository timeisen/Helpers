import sys

#Version 2 does not round the numbers.
alex_st_file=open(sys.argv[1],'r')
output_file=open(sys.argv[2],'w+')
for line in alex_st_file:
	gene=line.split("\t")[0]
	st_list=line.split("\t")[1].strip().split(",")
	if len(st_list)>=1:
		for tail in st_list[:-1]:
			print tail
			newline=gene+"\t"+str(float(tail))+"\n"
			sys.exit()
			output_file.write(newline)
alex_st_file.close()
output_file.close()
