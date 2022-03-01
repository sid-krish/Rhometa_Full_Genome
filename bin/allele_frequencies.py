#!/usr/bin/env python
import sys

import numpy as np
from Bio import SeqIO

inputFile = sys.argv[1]

# inputFile = "../Output/rho_15_sam_10_gen_20000/rho_15_sam_10_gen_20000_seqgenOut.fa"

allele_frequencies = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

records = list(SeqIO.parse(inputFile, "fasta"))

for record in records:
    seq = np.array(record.seq)

    unique, counts = np.unique(seq, return_counts=True)

    for unq, freq in zip(unique, counts):
        allele_frequencies[unq] = allele_frequencies[unq] + freq

with open('allele_frequencies.txt', 'w') as fileOut:
    dict_to_list = list(allele_frequencies.items())
    dict_to_list.sort()
    fileOut.write(f"Allele,Frequency\n")
    
    for allele, frequency in dict_to_list:
        fileOut.write(f"{allele},{frequency}")
        fileOut.write('\n')
    # convert later using dict(line)
