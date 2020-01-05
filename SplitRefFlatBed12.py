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
	handle.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(\
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

	txstart       = int(line.split("\t")[4])
	txEnd         = int(line.split("\t")[5])
	cdsStart      = int(line.split("\t")[6])
	cdsEnd        = int(line.split("\t")[7])
	exonStarts    = line.split("\t")[9].split(",")[:-1] #remove last element
	exonEnds      = line.strip().split("\t")[10].split(",")[:-1]

	if strand == ('-'): pos = False

	#Flags for parsing exons. 
	FpUtr = True
	Cds   = False
	TpUtr = False
	
	#Lists for bed12
	FpUtrLength, CdsLength, TpUtrLength = [], [], []
	FpUtrStart, CdsStart, TpUtrStart = [], [], []

	for exonIdx in range(len(exonStarts)): #the main code block for output exons
		#Tuple of starts and ends.
		currentExon = (int(exonStarts[exonIdx]), int(exonEnds[exonIdx]))

		#If the exon contains the cdsStart
		if FpUtr and currentExon[0] < cdsStart < currentExon[1] :
			FpUtrLength.append(cdsStart - currentExon[0])
			FpUtrStart.append(currentExon[0] - txstart)
			CdsLength.append(currentExon[1] - cdsStart)
			CdsStart.append(cdsStart - cdsStart)

			FpUtr = False
			Cds = True

		#Deals with an exon ending on a cdsStart. 
		elif cdsStart == currentExon[1]:
			FpUtrLength.append(cdsStart - currentExon[0])
			FpUtrStart.append(currentExon[0] - txstart)

			FpUtr = False
			Cds = True


		#If the exon is fully in the 5p utr
		elif FpUtr:
			FpUtrLength.append(currentExon[1] - currentExon[0])
			FpUtrStart.append(currentExon[0] - txstart)

		#If the exon is fully in the cds
		elif Cds and not currentExon[0] < cdsEnd < currentExon[1]:
			CdsLength.append(currentExon[1] - currentExon[0])
			CdsStart.append(currentExon[0] - cdsStart)

			#Deals with an exon ending on a cdsEnd. 
			if cdsEnd == currentExon[1]:
				Cds = False

		#If the exon contains the cdsEnd
		elif Cds and (currentExon[0] < cdsEnd < currentExon[1]):
			CdsLength.append(cdsEnd - currentExon[0])
			CdsStart.append(currentExon[0] - cdsStart)
			TpUtrLength.append(currentExon[1] - cdsEnd)
			TpUtrStart.append(cdsEnd - cdsEnd)
			
			Cds = False
			TpUtr = True

		#If the exon is fully in the 3p utr
		else:
			TpUtrLength.append(currentExon[1] - currentExon[0])
			TpUtrStart.append(currentExon[0] - cdsEnd)
		
	#writing files. 
	if pos:
		bed12writer(FpUtrBed, accession, strand, chrom, txstart, cdsStart,\
			FpUtrLength, FpUtrStart)
		bed12writer(CdsBed, accession, strand, chrom, cdsStart, cdsEnd,\
			CdsLength, CdsStart)
		bed12writer(TpUtrBed, accession, strand, chrom, cdsEnd, txEnd,\
			TpUtrLength, TpUtrStart)
	else: #just switches the 5p and 3p file writing. 
		bed12writer(TpUtrBed, accession, strand, chrom, txstart, cdsStart,\
			FpUtrLength, FpUtrStart)
		bed12writer(CdsBed, accession, strand, chrom, cdsStart, cdsEnd,\
			CdsLength, CdsStart)
		bed12writer(FpUtrBed, accession, strand, chrom, cdsEnd, txEnd,\
			TpUtrLength, TpUtrStart)


#Close the files
RefFlat.close()
FpUtrBed.close()
TpUtrBed.close()
CdsBed.close()