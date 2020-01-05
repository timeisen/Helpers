#!/bin/bash

#Takes as input a tab delimited file with sequence, then standard and greps the second argument for all lines in the first file. 
#Output is the third file. WARNING: Output file will be appended, not overwritten.

while read p; do
	std=$( echo "$p" | cut -f 2 )
	seq=$( echo "$p" | cut -f 1 )
	echo -e $std"\t"$( grep -c $seq $2 ) >>$3
done <$1