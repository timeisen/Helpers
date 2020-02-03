####################################
#This code will iterate over each codon in an orf and create a dictionary of codon frequencies for that orf.
#Need dictionary of dictionaries
#Timothy Eisen, 2015-4-24
####################################

from collections import defaultdict as dd
import sys, getopt
from Bio import SeqIO
from Bio.Data.CodonTable import unambiguous_dna_by_id

def main(argv): #accepts user inputs. Header function no yet defined
    try: opts, args=getopt.getopt(argv, "hi:o:c:", ["--help", "--header:"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    header=0
    for o, a in opts:
        if o=="-i": inputFile = a
        if o=="-o": outputFile = a
        if o=="--header": header=a
        if o=="-c": codon_text_file=a
        if o=="-h":
            usage()
    try:
        inputFile, outputFile, codon_text_file
    except: usage()
    (input, output, codon_list) = loadFiles(inputFile, outputFile, codon_text_file)
    codon_dictionary = dd(lambda: dd(float)) #change this to float before using frequencies
    for record in SeqIO.parse(inputFile, 'fasta'):
        codonExtraction(record, codon_dictionary)
    
    mRNA_codon_frequency(codon_dictionary, codon_list, output)
    
    close_files(input, output)

def loadFiles(inputFile, outputFile, codon_text_file):
    input=open(inputFile, 'rU')
    output=open(outputFile, 'w+')
    return(input, output)

def usage(): #if user inputs are not defined
    print '\nThis script calculates codon frequency for each gene in the FASTA file\n'
    print '-i\t\tInput file in FASTA format\n\n-o\t\tOutput file in tab delimited format \n\n-c\t\tCodon to calculate\n\n-h\t\tHelp\n\n--header\tTrims n lines of header from the input file'
    sys.exit(2)

def loadFiles(inputFile, outputFile, codon_text_file):
    input=open(inputFile, 'r')
    output=open(outputFile, 'w+')
    codons=open(codon_text_file, 'r')
    
    codon_list = []
    
    for line in codons: codon_list.append(line[8:11])
    codons.close()
    return(input, output, codon_list)

def codonExtraction(record, codon_dictionary):
    codon_counter = 0

    for i in range(0, len(record.seq), 3):
        codon = str(record.seq[i:i+3])
        codon_counter += 1
        codon_dictionary[record.id][codon]+=1
    for key, val in codon_dictionary[record.id].iteritems():
        codon_dictionary[record.id][key] = float(val)/float(codon_counter)

def mRNA_codon_frequency(codon_dictionary, codon_list, outfile):
	outfile.write('Gene_name\t%s\n'%('\t'.join(map(str,codon_list))))
	for gene, codons in codon_dictionary.iteritems():
		freq_list=[]
		freq_list.append(gene)
		for query_codon in codon_list:
			if query_codon in codons: freq_list.append (codons[query_codon])
			else: freq_list.append(0)
		outfile.write('%s\n'%('\t'.join(map(str,freq_list)))) 
		
def close_files(input, output):
	input.close()
	output.close()
		
		
if __name__ == "__main__":
    main(sys.argv[1:]) #note that this is only the arguments, not the script name (which is argv[0])
