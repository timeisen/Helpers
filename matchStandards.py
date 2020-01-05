#Needs to be run in a py3 virtual environment
#source /lab/solexa_bartel/teisen/RNAseq/virtualenv/py3/bin/activate

import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from collections import defaultdict as dd
import regex
import time

t0 = time.time()

stds = open(sys.argv[1],'r')
outFile = open(sys.argv[3],'w+')

stdDict = dd(str)
for line in stds:
	stdDict[line.split("\t")[1].strip()] = regex.compile('('+line.split("\t")[0]+'){e<=1}')


stdCount = dd(lambda: 0)
counter = 0
for record in SeqIO.parse(sys.argv[2], "fastq"):
    counter += 1
    if counter%100000==0:
    	t1 = time.time() 
    	print("%s reads processed in %s seconds per block ..." %(counter,round(t1-t0)) )
    	t0 = time.time()
    for std,stdseq in stdDict.items():
    	# if pairwise2.align.localms(record.seq,std,1,0,-100,-100,score_only=True)<19: continue
    	if stdseq.search(str(record.seq)):
    		stdCount[std] += 1
    		break

for k,v in stdCount.items():
	newline = k+"\t"+str(v)+"\n"
	outFile.write(newline)

stds.close()
outFile.close()