import fileinput
import sys

def reverse_complement(seq):
    """Get reverse complement of sequence"""
    nt_dict = {'A':'T', 'T': 'A', 'C': 'G', 'G':'C', 'N': 'N'}
    return ''.join([nt_dict[nt] for nt in seq][::-1])


for line in fileinput.input():
	newline = reverse_complement(line.strip()) + "\n"
	sys.stdout.write(newline)
