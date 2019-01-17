
import sys
import os
from Bio import SeqIO
from collections import Counter
import pandas as pd
import numpy as np
from multiprocessing import Pool

dc_nt = {'A': 1, 'T': 2, 'C': 3, 'G': 4}

def nt2int(nt,missing_data = -9):
    '''
    convert nt to int based on dc_nt. missing_data will be -0
    '''
    if nt in dc_nt:
        return dc_nt[nt]
    return missing_data

def fasta2dataframeOneThread(e, missing_data = -9):
    seq = str(e.seq)
    seqInt = np.array([nt2int(e,missing_data) for e in seq], dtype = np.int8)
    return seqInt

def fasta2dataframe(filename, missing_data = -9, threads = 32):
    '''
    filename is a fasta file
    return a dataframe, with seq.id as index
    '''
    ls_fasta = list(SeqIO.parse(filename,'fasta'))
    pool = Pool(threads)
    ls_int = pool.starmap(fasta2dataframeOneThread, [(e,missing_data) for e in ls_fasta])
    pool.close()
    df_nt = pd.DataFrame()
    for n, e in enumerate(ls_fasta):
        seqInt = ls_int[n]
        df_nt[e.id] = seqInt
    return df_nt

def multipleFastaHaploid2STRUCTUREinputOneThread(df_nt, gapcut = 0.8, Ncut=4, missing_data = -9):
    '''
    df_nt is a dataframe of nt in np.int8. 
    return a dataframe based on gapcut = 0.8, Ncut=4, missing_data = -9
    '''
    df_nt = df_nt.copy() # make a copy to avoid messing up
    seq_len, seq_count = df_nt.shape
    max_missing = seq_count * (1-gapcut)
    rows_keep = []
    rows_index_keep = []
    # change based on Ncut. Any with count less than Ncut will be changed to -9, missing_data
    # keep sites based on gapcut and base types should be greater than 1
    for row, row_value in df_nt.iterrows():
        dc_row = dict(row_value)
        ntCounts = Counter(dc_row.values())
        for col in dc_row:
            nt = dc_row[col]
            if ntCounts[nt] < Ncut and nt != -9:
                dc_row[col] = missing_data
        ntCounts = Counter(dc_row.values())
        if missing_data in ntCounts:#if -9 in ntCounts, than at least three different keys in ntCounts
            if ntCounts[missing_data] < max_missing and len(ntCounts)>=3:
                rows_keep.append(dc_row)
                rows_index_keep.append(row)
        else:# -9 not in ntCounts, at least two kinds of keys in ntCounts
            if len(ntCounts) >= 2:
                rows_keep.append(dc_row)
                rows_index_keep.append(row)
    
    df_snp = pd.DataFrame(rows_keep,index=rows_index_keep)
    return df_snp

def multipleFastaHaploid2STRUCTUREinput(filename,outfile = 'default', gapcut = 0.8, Ncut=4, missing_data = -9, threads = 32, header_map_dist=True):
    '''
    filename is a fasta file with haploid sequences.
    output a file for STRUCTURE to use as input
    gapcut is the ratio of allowed gaps in each position
    Ncut is the minimum required counts for each element at each positon to shown in the output
    work for nucleotide sequences
    if outfile is None, return the result as dataframe
    header_map_dist is the header line indicating Map Distances between Loci. If True, include header_map_dist in writing output. Else, do not include
    '''
    
    # store the input fasta file to a dataframe. Convert 'ATCG' to numbers 1234
    df_nt = fasta2dataframe(filename,missing_data=missing_data, threads=threads)
    seq_len, seq_count = df_nt.shape
    print('finished reading fasta file. seq_len, seq_count are ',df_nt.shape)
    
    df_nt_splits = [df_nt.loc[e] for e in np.array_split(df_nt.index,threads)]
    pool = Pool(threads)
    df_snp_splits = pool.starmap(multipleFastaHaploid2STRUCTUREinputOneThread, [(e,gapcut, Ncut, missing_data) for e in df_nt_splits])
    pool.close()
    df_snp = pd.concat(df_snp_splits)
    
    df_snp = df_snp.T
    print('sequence count', df_snp.shape[0], 'SNP count', df_snp.shape[1]+1)
    
    df_snp[-1]=0
    df_snp=df_snp.sort_index(axis=1)
    
    header_col = [-1]
    header_col.append(df_snp.columns[1])
    for i in range(1,len(df_snp.columns)-1):
        header_col.append(df_snp.columns[i+1] - df_snp.columns[i])
    df_snp.columns = header_col
    
    if outfile == 'False':
        return df_snp
    if outfile == 'default':
        outfile = filename + '.STRUCTURE'
    
    if header_map_dist:
        df_snp.to_csv(outfile,sep='\t')
    else:
        df_snp.to_csv(outfile,sep='\t',header=False)
    print('finish wrinting file', outfile)
    print('sequence count', df_snp.shape[0], 'SNP count', df_snp.shape[1])

description = '''
    filename is a fasta file with haploid sequences.
    output a file for STRUCTURE to use as input
    gapcut is the ratio of allowed gaps in each position
    Ncut is the minimum required counts for each element at each positon to shown in the output
    work for nucleotide sequences
    if outfile is None, return the result as dataframe
    default thread is 32 for multiple processing
    header_map_dist is the header line indicating Map Distances between Loci. If True, include header_map_dist in writing output. Else, do not include

'''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file fasta file', required=True)
    parser.add_argument('-o','--output',help = 'output file to store the STRUCTURE input file. default, input filename.STRCUTRUE', default='default')
    parser.add_argument('-t','--threads',help = 'number of threads to use', default = 32, type=int)
    parser.add_argument('-g','--gapcut',help = 'gapcut ratio, minimum ratio of bases in each site, default 0.8', default = 0.8,type=float)
    parser.add_argument('-N','--Ncut',help = 'minumum number of counts to consider as a SNP, default 4. (>=4)',default = 4, type=int)
    parser.add_argument('-M','--missing_data',help = 'set value other than ATCG to this value, default -9',default = -9, type=int)
    parser.add_argument('-H','--header_map_dist',help = 'whether include header_map_dist in output file, 1 for True, 0 for False. default 1',default = 1,type=bool, choices=[0,1])
    f = parser.parse_args()
    multipleFastaHaploid2STRUCTUREinput(filename = f.input,outfile = f.output, gapcut = f.gapcut, Ncut=f.Ncut, missing_data = f.missing_data, threads=f.threads,header_map_dist=f.header_map_dist)
    




