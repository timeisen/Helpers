#this script will take as input a FASTA file and output the length of all the sequences
#the third argument is a boolean specifying whether the header is included

import Bio, sys

from Bio import SeqIO

transcripts = open(sys.argv[1],"rU")
outputfile = open(sys.argv[2],"w+")
if sys.argv[3] == "True": outputfile.write("accession\tlength\n") #add a header

for record in SeqIO.parse(transcripts, "fasta"):
    outputfile.write("%s\t%s\n"%(record.id,len(record.seq)))

outputfile.close()