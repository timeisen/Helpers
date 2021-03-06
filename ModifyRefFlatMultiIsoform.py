################################################################################
# ModifyRefFlatMultiIsoform.py
# Timothy J Eisen
# 2020 02 03
# Takes in a refFlat file and a list of 3' end annotations from PAL-seq data
# Modifies the txEnd coordinates of the refFlat file, keeping genes that are not
# in the PAL-seq data unchanged. 
# Outputs multiple entries if they exist in the PAL-seq annotations. 
# Format of input annotation: accession\tchrom\tstrand\t3pEndPos\ttagFraction\n
################################################################################


import sys, pdb
from collections import defaultdict as dd

RefFlat  = open(sys.argv[1],'r')
PalAnno = open(sys.argv[2],'r')
RefFlatOut   = open(sys.argv[3],'w+')

#Skip the header
next(PalAnno)

#Iterate through the pal-seq annotations, creating a dictionary of all 3p ends
# for each gene
AnnoDict = dd(lambda: ([],[]))
for line in PalAnno:
	#Parse the line
	accession       = line.split("\t")[0]
	chrom		    = line.split("\t")[1]
	strand          = line.split("\t")[2]
	txEnd           = int(line.split("\t")[3])
	tagFrac         = float(line.split("\t")[4].strip())

	AnnoDict[accession][0].append(txEnd)
	AnnoDict[accession][1].append(tagFrac)

#iterate through RefFlat
for line in RefFlat:


	#Parse the line.
	symbol        = line.split("\t")[0]
	strand        = line.split("\t")[3]
	chrom         = line.split("\t")[2]
	accession     = line.split("\t")[1]
	NumberOfExons = int(line.split("\t")[8])
	txStart       = int(line.split("\t")[4])
	txEnd         = int(line.split("\t")[5])
	cdsStart      = int(line.split("\t")[6])
	cdsEnd        = int(line.split("\t")[7])
	exonStarts    = line.split("\t")[9].split(",")[:-1] #remove last element
	exonEnds      = line.strip().split("\t")[10].split(",")[:-1]

	if accession not in AnnoDict:
		ListOfLine = line.split("\t")
		ListOfLine[1] = ListOfLine[1] +":0:1"
		RefFlatOut.write("\t".join(ListOfLine))
		continue

	for TxAnnoIdx, TxEndAnno in enumerate(AnnoDict[accession][0]):
		accession_mod = accession + ":" + str(TxAnnoIdx) + ":" + \
			str(AnnoDict[accession][1][TxAnnoIdx])
		txEndMod = TxEndAnno
		exonStartsMod = exonStarts[:]
		exonEndsMod = exonEnds[:]
	
		if strand == '+':
			#only for positive right now.
			for exonIdx in range(len(exonStarts) - 1,0,-1): #the main code block for output exons
				#Tuple of starts and ends.
				currentExon = (int(exonStarts[exonIdx]), int(exonEnds[exonIdx]))
				if currentExon[0] < txEndMod: #or =?
					exonEndsMod[exonIdx] = str(txEndMod)
					ListEnd = exonIdx
					break
				# else: 
				# 	del exonStartsMod[exonIdx]
				# 	del exonEndsMod[exonIdx]
			exonStartsMod = exonStartsMod[:(ListEnd + 1)]
			exonEndsMod = exonEndsMod[:(ListEnd + 1)]
		
			if cdsEnd > txEndMod:
				cdsEnd = txEndMod
				#This warning is all I have for issues like this, but it could be expanded. 
				print("Warning: {} contains no 3p UTR in modified annotations."\
					.format(accession_mod))
			txEnd = txEndMod
	
	
		elif strand == '-':
			for exonIdx in range(len(exonStarts)): #the main code block for output exons
				#Tuple of starts and ends.
				currentExon = (int(exonStarts[exonIdx]), int(exonEnds[exonIdx]))
				if txEndMod < currentExon[1]: #or =?
					exonStartsMod[exonIdx] = str(txEndMod)
					ListStart = exonIdx
					break
				# else: 
				# 	del exonStartsMod[exonIdx]
				# 	del exonEndsMod[exonIdx]
		
			exonStartsMod = exonStartsMod[ListStart:]
			exonEndsMod = exonEndsMod[ListStart:]
	
			if cdsStart < txEndMod:
				cdsStart = txEndMod
				print("Warning: {} contains no 3p UTR in modified annotations."\
					.format(accession_mod))
			txStart = txEndMod
	
		NewRefFlatLine = "\t".join([symbol, accession_mod, chrom, strand, \
			str(txStart), str(txEnd), str(cdsStart), str(cdsEnd), \
			str(len(exonStartsMod)), \
			",".join(exonStartsMod) + ",", ",".join(exonEndsMod) + ","])
		RefFlatOut.write(NewRefFlatLine + "\n")


#Close the files
RefFlat.close()
PalAnno.close()
RefFlatOut.close()
