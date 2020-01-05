#position_coordinate

import sys
from collections import defaultdict as dd

input_reads=open(sys.argv[1],'r') #of the form 1101:3453:23424 for example
output_file=open(sys.argv[2],'w+')

output_file.write("x\ty\tz\tcount\n")

position_coordinate=dd(lambda: 0)
for line in input_reads: 
	z=int(list(line.split(":")[0])[0])
	x=int(list(line.split(":")[0])[1])
	y=int("".join(list(line.split(":")[0])[2:4]))
	key=(x,y,z)

	position_coordinate[key]+=1

for  key, val in position_coordinate.iteritems():
	outline='\t'.join('%s' % x for x in key)+"\t"+str(val)+"\n"
	output_file.write(outline)

input_reads.close()
output_file.close()


##__________________________R Commands__________________________#

#library(ggplot2)
#p1<-ggplot(df,aes(x,y))+geom_raster(aes(fill=count),interpolate=TRUE)

