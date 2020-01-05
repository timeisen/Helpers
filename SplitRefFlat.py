################################################################################
# SplitRefFlat.py
# Timothy J Eisen
# 2020 01 04
# Creates 3 bed files from a refflat file consisting of 3' UTR, CDS, and 5' UTR
#     sequences. These can be used to get the corresponding fasta data.
# Takes into account strandedness and exon junctions.
# The resulting bed files are ordered. 
################################################################################

#Could I have just ignored strand information and swap the UTR output files for
# negative stranded genes?	
#Maybe, but the bed file order would be wrong. It might be nice to preserve this.

import sys, pdb

RefFlat  = open(sys.argv[1],'r')
FpUtrBed = open(sys.argv[2],'w+')
CdsBed   = open(sys.argv[3],'w+')
TpUtrBed = open(sys.argv[4],'w+')

for line in RefFlat:
	pos = True

	#Parse the line.
	strand        = line.split("\t")[3]
	chrom         = line.split("\t")[2]
	accession     = line.split("\t")[1]
	NumberOfExons = int(line.split("\t")[8])

	if strand == "+":
		txstart       = int(line.split("\t")[4])
		txEnd         = int(line.split("\t")[5])
		cdsStart      = int(line.split("\t")[6])
		cdsEnd        = int(line.split("\t")[7])
		exonStarts    = line.split("\t")[9].split(",")[:-1] #remove last element
		exonEnds      = line.strip().split("\t")[10].split(",")[:-1]

	elif strand == "-": #Swap the starts and ends if the gene is - stranded
		pos           = False
		txstart       = int(line.split("\t")[5])
		txEnd         = int(line.split("\t")[4])
		cdsStart      = int(line.split("\t")[7])
		cdsEnd        = int(line.split("\t")[6])
		exonStarts    = line.strip().split("\t")[10].split(",")[-2::-1]
		exonEnds      = line.split("\t")[9].split(",")[-2::-1]

	else:
		raise ValueError('Strand parsing failure.')

	#Flags for parsing exons. 
	FpUtr = True
	Cds   = False
	TpUtr = False
	for exonIdx in range(len(exonStarts)): #the main code block for output exons
		#Tuple of starts and ends.
		currentExon = (int(exonStarts[exonIdx]), int(exonEnds[exonIdx]))

		#If the exon contains the cdsStart
		if (currentExon[0] < cdsStart < currentExon[1] and pos) or \
			(currentExon[0] > cdsStart > currentExon[1] and not pos) and FpUtr:
			if pos:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx],cdsStart, accession, strand))
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, cdsStart,exonEnds[exonIdx], accession, strand))
			else:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, cdsStart, exonStarts[exonIdx], accession, strand))
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx], cdsStart, accession, strand))

			FpUtr = False
			Cds = True

		#Deals with an exon ending on a cdsStart. 
		elif (cdsStart == currentExon[0] and not pos) or\
			(cdsStart == currentExon[1] and pos):
			if pos:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx],\
					exonEnds[exonIdx], accession, strand))
			else:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx],\
					exonStarts[exonIdx], accession, strand))
			
			FpUtr = False
			Cds = True


		#If the exon is fully in the 5p utr
		elif FpUtr:
			if pos:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx],\
					exonEnds[exonIdx], accession, strand))
			else:
				FpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx],\
					exonStarts[exonIdx], accession, strand))

		#If the exon is fully in the cds
		elif Cds and not (currentExon[0] < cdsEnd < currentExon[1] and pos) \
			and not (currentExon[0] > cdsEnd > currentExon[1] and not pos):
			if pos:
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx],\
					exonEnds[exonIdx], accession, strand))
			else:
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx],\
					exonStarts[exonIdx], accession, strand))
			#Deals with an exon ending on a cdsEnd. 
			if (cdsEnd == currentExon[0] and not pos) or\
				(cdsEnd == currentExon[1] and pos):
				Cds = False

		#If the exon contains the cdsEnd
		elif Cds and (currentExon[0] < cdsEnd < currentExon[1] and pos) or \
			(currentExon[0] > cdsEnd > currentExon[1] and not pos):
			if pos:
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx], cdsEnd, accession, strand))
				TpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, cdsEnd,exonEnds[exonIdx], accession, strand))
			else: 
				CdsBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, cdsEnd, exonStarts[exonIdx], accession, strand))
				TpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx], cdsEnd, accession, strand))

			Cds = False
			TpUtr = True

		#If the exon is fully in the 3p utr
		else:
			if pos:
				TpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonStarts[exonIdx],\
					exonEnds[exonIdx], accession, strand))
			else:
				TpUtrBed.write("{}\t{}\t{}\t{}\t0\t{}\n".format(\
					chrom, exonEnds[exonIdx],\
					exonStarts[exonIdx], accession, strand))				




#Close the files
RefFlat.close()
FpUtrBed.close()
TpUtrBed.close()
CdsBed.close()