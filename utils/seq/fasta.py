# -*- coding: utf-8 -*-
"""
functions to process fasta file, fasta sequences
Created on Mon Dec  3 17:19:22 2018

@author: ATPs
"""

from collections import Counter
from Bio import SeqIO
import pandas as pd

def _baseCout(seq):
    '''
    Helper function
    count the number of bases in a seq. seq is string
    return a dictionary, with bases and their counts
    '''
    return dict(Counter(seq))

def baseCout(seqs, outfile=None):
    '''
    seqs is a filename, a list of sequences, or a dictionary of sequences.
    if seqs is dictionary, use keys of dictionary as sequence names
    if seqs is a list:
        if sequences is SeqIO.SeqRecord, use id of sequences as names
        if sequences is str, use 0,1,2,... as id of sequences
    return a dataframe of basecounts in each seqs
    if outfile is not None, write a csv file of returned dataframe to outfile
    '''
    # get sequence names and sequences in str
    if isinstance(seqs, dict):
        print('input seqs is a dict')
        names = seqs.keys()
        sequences = seqs.values()
    elif isinstance(seqs,list):
        if isinstance(seqs[0], str):
            print('input seqs is a list of str')
            names = list(range(len(seqs)))
            sequences = seqs
        elif isinstance(seqs[0], SeqIO.SeqRecord):
            print('input seqs is a list of SeqIO.SeqRecord')
            names = [e.id for e in seqs]
            sequences = [str(e.seq) for e in seqs]
        else:
            print('input seqs is a list, but not in the right format')
            return None
    elif isinstance(seqs, str):
        print('input seqs is a filename of sequences in fasta format')
        seqs = list(SeqIO.parse(seqs,'fasta'))
        names = [e.id for e in seqs]
        sequences = [str(e.seq) for e in seqs]
    else:
        print('wrong input format for seqs')
    
    # Count bases for seqs
    results = []
    for name, sequence in zip(names, sequences):
        dc_bases = _baseCout(sequence)
        dc_bases['name'] = name
        results.append(dc_bases)
    
    # return dataframe
    df = pd.DataFrame(results)
    df = df.set_index('name')
    if outfile is not None:
        df.to_csv(df)
    return df

def alignmentFilterGapSeqs(seqs_in, seqs_out=None, max_gapratio = 0.8):
    '''
    seqs_in is a list of SeqIO.SeqRecord, or a filename of fasta file
    if seqs_out is not None, write the result to seqs_out
    return a list of sequences in SeqIO.SeqRecord format
    '''
    if isinstance(seqs_in,str):
        seqs_in = list(SeqIO.parse(seqs_in,'fasta'))
    
    ls_out = []
    for seq in seqs_in:
        seqlen = len(seq.seq)
        max_gap_count = max_gapratio * seqlen
        dc_bases = _baseCout(str(seq.seq))
        if dc_bases['-'] < max_gap_count:
            ls_out.append(seq)
    
    if seqs_out is not None:
        fout = open(seqs_out,'w')
        for seq in ls_out:
            fout.write('>'+seq.id+'\n'+str(seq.seq)+'\n')
        fout.close()
    
    return ls_out

