# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 14:39:41 2018

@author: ATPs
"""

import sys
filename = sys.argv[1]
from Bio import SeqIO
fout = open('genome.scf_size','w')

for s in SeqIO.parse(filename,'fasta'):
    seq = str(s.seq)
    seq.upper()
    seqlen = len(seq)
    seqlenNoN = seqlen - seq.count('N')
    fout.write('{}\t{}\t{}\n'.format(s.id, seqlen, seqlenNoN))

fout.close()