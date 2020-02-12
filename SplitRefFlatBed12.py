################################################################################
# SplitRefFlat.py
# Timothy J Eisen
# 2020 01 05
# Creates 3 bed12 files from a refflat file consisting of 3' UTR, CDS, and 5' UTR
#     sequences. These can be used to get the corresponding fasta data.
# Takes into account strandedness and exon junctions.
# 
################################################################################


import sys, pdb

RefFlat  = open(sys.argv[1],'r')
FpUtrBed = open(sys.argv[2],'w+')
CdsBed   = open(sys.argv[3],'w+')
TpUtrBed = open(sys.argv[4],'w+')

def bed12writer(handle, accession, strand, chrom, start, \
	end, featSizes,featStarts):
	if len(featSizes) == 0: pass #Don't write lines that don't have features. 
	else: handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
		chrom, start, end, accession, 0, strand, start, end, "255,0,0", \
		len(featSizes), ",".join([str(i) for i in featSizes]),\
		",".join([str(i) for i in featStarts])))

for line in RefFlat:

	pos = True

	#Parse the line.
	strand        = line.split("\t")[3]
	chrom         = line.split("\t")[2]
	accession     = line.split("\t")[1]
	NumberOfExons = int(line.split("\t")[8])

	TxStart       = int(line.split("\t")[4])
	TxEnd         = int(line.split("\t")[5])
	CdsStart      = int(line.split("\t")[6])
	CdsEnd        = int(line.split("\t")[7])
	ExonStarts    = line.split("\t")[9].split(",")[:-1] #remove last element
	ExonEnds      = line.strip().split("\t")[10].split(",")[:-1]

	if strand == ('-'): pos = False

	#Flags for parsing exons. 
	FpUtr = True
	Cds   = False
	TpUtr = False
	
	#Lists for bed12
	FpUtrLengthList, CdsLengthList, TpUtrLengthList = [], [], []
	FpUtrStartList, CdsStartList, TpUtrStartList = [], [], []

	#For dealing with a 3' UTR being its own exon.
	TpUtrBegin = CdsEnd
	TpUtrNotContig = False

	for exonIdx in range(len(ExonStarts)): #the main code block for output exons
		#Tuple of starts and ends.
		currentExon = (int(ExonStarts[exonIdx]), int(ExonEnds[exonIdx]))

		#If the exon contains the CdsStart
		if FpUtr and currentExon[0] < CdsStart < currentExon[1] :
			FpUtrLengthList.append(CdsStart - currentExon[0])
			FpUtrStartList.append(currentExon[0] - TxStart)
			CdsLengthList.append(currentExon[1] - CdsStart)
			CdsStartList.append(CdsStart - CdsStart)

			FpUtr = False
			Cds = True

		#Deals with CdsStart and TxStart being the same. (no 5p UTR)
		elif CdsStart == currentExon[0]:
			CdsLengthList.append(currentExon[1] - CdsStart)
			CdsStartList.append(CdsStart - CdsStart)

			FpUtr = False
			Cds = True

		#Deals with an exon ending on a CdsStart. 
		elif CdsStart == currentExon[1]:
			FpUtrLengthList.append(CdsStart - currentExon[0])
			FpUtrStartList.append(currentExon[0] - TxStart)

			FpUtr = False
			Cds = True

		#If the exon is fully in the 5p utr
		elif FpUtr:
			FpUtrLengthList.append(currentExon[1] - currentExon[0])
			FpUtrStartList.append(currentExon[0] - TxStart)

		#If the exon is fully in the cds
		elif Cds and not (currentExon[0] < CdsEnd <= currentExon[1]):
			CdsLengthList.append(currentExon[1] - currentExon[0])
			CdsStartList.append(currentExon[0] - CdsStart)

		#If the exon contains the CdsEnd
		elif Cds and (currentExon[0] < CdsEnd <= currentExon[1]):
			CdsLengthList.append(CdsEnd - currentExon[0])
			CdsStartList.append(currentExon[0] - CdsStart)
			Cds = False
			TpUtr = True

			#Deals with a 3' UTR being its own exon. 
			if CdsEnd == currentExon[1]: 
				TpUtrNotContig = True 

			#Deals with an exon ending on a CdsEnd. 
			elif CdsEnd == TxEnd: TpUtr = False
			else: 
				TpUtrLengthList.append(currentExon[1] - CdsEnd)
				TpUtrStartList.append(CdsEnd - CdsEnd)
			


		#If the exon is fully in the 3p utr
		elif TpUtr:
			if TpUtrNotContig: 
				TpUtrBegin = currentExon[0]
				TpUtrNotContig = False
			# 	continue
			# 	# pdb.set_trace()
			TpUtrLengthList.append(currentExon[1] - currentExon[0])
			TpUtrStartList.append(currentExon[0] - TpUtrBegin)

		
	#writing files.  
	if pos:
		bed12writer(FpUtrBed, accession, strand, chrom, TxStart, CdsStart,\
			FpUtrLengthList, FpUtrStartList)
		bed12writer(CdsBed, accession, strand, chrom, CdsStart, TpUtrBegin,\
			CdsLengthList, CdsStartList)
		bed12writer(TpUtrBed, accession, strand, chrom, TpUtrBegin, TxEnd,\
			TpUtrLengthList, TpUtrStartList)
	else: #just switches the 5p and 3p file writing. 
		bed12writer(TpUtrBed, accession, strand, chrom, TxStart, CdsStart,\
			FpUtrLengthList, FpUtrStartList)
		bed12writer(CdsBed, accession, strand, chrom, CdsStart, TpUtrBegin,\
			CdsLengthList, CdsStartList)
		bed12writer(FpUtrBed, accession, strand, chrom, TpUtrBegin, TxEnd,\
			TpUtrLengthList, TpUtrStartList)


#Close the files
RefFlat.close()
FpUtrBed.close()
TpUtrBed.close()
CdsBed.close()