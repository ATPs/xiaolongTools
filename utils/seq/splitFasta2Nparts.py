# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 11:14:48 2018

@author: ATPs
"""
from Bio import SeqIO
import os

def splitFasta2Nparts(seqs, N=10, output=None, method = 1):
    '''
    seqs is a filename, or a list of sequences
    if output is None
        if seqs is a filename, ouput is filename+'.splitN'
        if seqs is a list, output is 'seqs.splitN'
    N is the number of files to return
    method to split the sequences. three choices, 0 for Continue, 1 for EqualStep, 2 for RandomChoice
    '''
    if isinstance(seqs,str):
        print('input is a filename')
        if output is None:
            output = seqs + '.split{}.'.format(N)
        seqs = SeqIO.parse(seqs,'fasta')
    elif isinstance(seqs,list):
        print('input is a list')
        if output is None:
            output = 'seqs.split{}.'.format(N)
    else:
        print('unknown input format for seqs')
    
    dirname = os.path.dirname(output)
    if dirname != '':
        if not os.path.exists(dirname):
            print('make folder', dirname)
            os.makedirs(dirname)
    
    ls_fout = [open(output+str(i), 'w') for i in range(N)]
    
    if method == 1:
        for n, seq in enumerate(seqs):
            ls_fout[n%N].write('>'+seq.description+'\n'+str(seq.seq)+'\n')
        for fout in ls_fout:
            fout.close()
        print('done!')
        return None
    print('method other the 1 is not working now')
    return None

description = '''
seqs is a filename, or a list of sequences
    if output is None
        if seqs is a filename, ouput is filename+'.splitN'
        if seqs is a list, output is 'seqs.splitN'
    N is the number of files to return
    method to split the sequences. three choices, 0 for Continue, 1 for EqualStep, 2 for RandomChoice
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-s','--seqs', help = 'seqs is a filename, or a list of sequences', required=True)
    parser.add_argument('-o','--output',help = 'where to store output file', required = False, default = None)
    parser.add_argument('-N','--number', help = 'number of output files. split to N parts', default = 10, type=int)
    parser.add_argument('-m','--method', help = 'method to split the sequences. three choices, 0 for Continue, 1 for EqualStep, 2 for RandomChoice', default = 1, type=int, choices = [0,1,2])
    f = parser.parse_args()
    splitFasta2Nparts(seqs=f.seqs, N=f.number, output=f.output, method = f.method)