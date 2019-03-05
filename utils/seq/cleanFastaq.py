# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 12:26:20 2019

@author: ATPs
"""


threads = 32
file_kmer1 = '/home/xcao/w/2018Calephelis/20190207cleanup15109F05/Arumecla_aruma.21mer.fa'
file_kmer2 = '/home/xcao/w/2018Calephelis/20190207cleanup15109F05/Calephelis_freemani.21mer.fa'
file_fq1 = '/home/xcao/w/2018Calephelis/20190207cleanup15109F05/15109F05_R1.fastq'
file_fq2 = '/home/xcao/w/2018Calephelis/20190207cleanup15109F05/15109F05_R2.fastq'
file_kmerT = '/home/xcao/w/2018Calephelis/20190207cleanup15109F05/target.21mer.fa'

from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import os
import glob
from multiprocessing import Pool
import pickle
import pandas as pd
import sqlite3

dcBase = {'A':'0','T':'1','C':'2', 'G':'3'}

def kmer2int(kmer):
    '''
    kmer is like 'ATAAAGC', convert it to int.
    A=0 T=1 C=2 G=3
    '''
    kmer = kmer.upper()
    l = [dcBase[e] for e in kmer]
    s = ''.join(l)
    return int(s,4)
    
    
def kmerfile2npArr(filename):
    '''
    given a file of fasta format of kmers, return two numpy array of int. one is change kmer to int, one is the count of each kmer.
    '''
    txt = open(filename).read()
    n = txt.count('>')
    del txt
    npKmer = np.zeros(n,dtype=np.int64)
    npCount = np.zeros(n,dtype=np.int32)
    for n, s in enumerate(SeqIO.parse(filename,'fasta')):
        npKmer[n] = kmer2int(str(s.seq))
        npCount[n] = int(s.id)
    return npKmer, npCount

def kmerfile2npArrThread(file_kmerfa, threads, outfile=False,tempFolderInMem=False):
    '''
    kmerfile2npArr work in multiple threads. Return a dataframe
    if outfile is True, save the files in binary format
    tempFolderInMem is whether to store the temporary files in memory. default False, the same location as input files
    '''
    if not tempFolderInMem:
        temp_kmerfa = file_kmerfa
    else:
        basename = os.path.basename(file_kmerfa)
        temp_kmerfa = '/dev/shm/'+basename
    
    cmd = 'split {file_kmerfa} -l 10000000 -a 5 -d {temp_kmerfa}_'.format(file_kmerfa=file_kmerfa,temp_kmerfa=temp_kmerfa)
    os.system(cmd)
    files = glob.glob('{temp_kmerfa}_*'.format(temp_kmerfa=temp_kmerfa))
    
    pool = Pool(threads)
    results = pool.map(kmerfile2npArr, files)
    pool.close()
    os.system('rm {temp_kmerfa}_*'.format(temp_kmerfa=temp_kmerfa))
    
    npKmer = np.concatenate([e[0] for e in results])
    npCount = np.concatenate([e[1] for e in results])
    
    df = pd.DataFrame(data={'kmer':npKmer,'count':npCount})
    del npKmer, npCount
    df = df.sort_values(by='kmer', kind='mergesort')
    df = df.reset_index(drop=True)
    npKmer = df['kmer'].values
    npCount = df['count'].values
    del df
    if outfile:
        fout = open(file_kmerfa+'.binary.kmer','wb')
        pickle.dump(npKmer,fout,protocol=4)
        fout.close()
        fout = open(file_kmerfa+'.binary.count','wb')
        pickle.dump(npCount,fout,protocol=4)
        fout.close()
    
    return npKmer, npCount

def seq2kmer(seq, kmerlen=21, bothstrand=True):
    '''
    seq is a sequence of ATCG, return a list of kmer with length of kmerlen
    if bothstrand, also return kmers from the reverse completement sequence
    '''
    seqlen = len(seq)
    if seqlen < kmerlen:
        return []
    kmers = []
    for i in range(0,seqlen-kmerlen+1):
        kmers.append(seq[i:i+kmerlen])
    if bothstrand:
        seq = Seq(seq)
        seq_rc = seq.reverse_complement()
        seq_rc = str(seq_rc)
        for i in range(0,seqlen-kmerlen+1):
            kmers.append(seq_rc[i:i+kmerlen])
    return kmers


def findIdenticalSitesOfQueryInSubject(query, subject):
    '''
    query, subject are numpy array
    return sites in subject that is identical to query
    '''
    locs = np.searchsorted(subject, query)
    query = query[locs < len(subject)]
    locs = locs[locs < len(subject)]
    return locs[query == subject[locs]]

def splitFq_each(file_fq1, file_fq2, file_using):
    '''
    split file_fq1, file_fq2 to fout1 and fout2 based on kmer info
    return a numpy array with 0 or 1, 2 to indicate which sample the reads belongs to. 
    0 means to remove
    1: npKmer1
    2: npKmer2
    '''
    npKmer_1, npCount_1, npKmer_2, npCount_2 = pickle.load(open(file_using,'rb'))
    seqs_1 = list(SeqIO.parse(file_fq1,'fastq'))
    seqs_2 = list(SeqIO.parse(file_fq2,'fastq'))
    
    seqCount_1 = len(seqs_1)
    seqCount_2 = len(seqs_2)
    if seqCount_1 != seqCount_2:
        print(file_fq1, file_fq2, 'not the same length')
        return None
    seqCount = len(seqs_1)
    result = np.zeros(seqCount, dtype=np.int8)
    
    N = -1
    for seq1, seq2 in zip(seqs_1, seqs_2):
        N += 1
        seq1seq = str(seq1.seq)
        seq2seq = str(seq2.seq)
        if 'N' in seq1seq or 'N' in seq2seq:
            continue
        kmers = seq2kmer(seq1seq) + seq2kmer(seq2seq)
        if len(kmers) == 0:
            continue
        kmerInt = np.array([kmer2int(e) for e in kmers], dtype=int)
        locs1 = findIdenticalSitesOfQueryInSubject(kmerInt, npKmer_1)
        locs2 = findIdenticalSitesOfQueryInSubject(kmerInt, npKmer_2)
        if len(locs1) == 0 and len(locs2) == 0:
            continue
#        score_1 = npCount_1[locs1].sum()
#        score_2 = npCount_2[locs2].sum()
        score_1 = len(locs1)
        score_2 = len(locs2)
        if score_1 >= score_2:
            result[N] = 1
        else:
            result[N] = 2
    
    #write files
    fout0_1 = open(file_fq1+'.0','w')
    fout0_2 = open(file_fq2+'.0','w')
    fout1_1 = open(file_fq1+'.1','w')
    fout1_2 = open(file_fq2+'.1','w')
    fout2_1 = open(file_fq1+'.2','w')
    fout2_2 = open(file_fq2+'.2','w')
    dc_out = {0:[fout0_1, fout0_2], 1:[fout1_1, fout1_2], 2:[fout2_1, fout2_2]}
    for seq1, seq2, n in zip(seqs_1, seqs_2, result):
        dc_out[n][0].write(seq1.format('fastq'))
        dc_out[n][1].write(seq2.format('fastq'))
    for v in dc_out.values():
        v[0].close()
        v[1].close()
    
    return result


def splitFq(file_fq1, file_fq2, file_kmer1, file_kmer2, file_kmerT, threads = 32):
    '''
    file_fq1, file_fq2 is pair end file
    file_kmer1, file_kmer2 is two kmer file
    '''
    # kmercount
    try:
        npKmer_1 = pickle.load(open(file_kmer1+'.binary.kmer', 'rb'))
        npCount_1 = pickle.load(open(file_kmer1+'.binary.count', 'rb'))
    except:
        print(file_kmer1, 'binary file not exist, create new ones')
        npKmer_1, npCount_1 = kmerfile2npArrThread(file_kmerfa=file_kmer1, threads=threads, outfile=True,tempFolderInMem=True)
    
    try:
        npKmer_2 = pickle.load(open(file_kmer2+'.binary.kmer', 'rb'))
        npCount_2 = pickle.load(open(file_kmer2+'.binary.count', 'rb'))
    except:
        print(file_kmer2, 'binary file not exist, create new ones')
        npKmer_2, npCount_2 = kmerfile2npArrThread(file_kmerfa=file_kmer2, threads=threads, outfile=True,tempFolderInMem=True)
    
    try:
        npKmer_t = pickle.load(open(file_kmerT+'.binary.kmer', 'rb'))
        npCount_t = pickle.load(open(file_kmerT+'.binary.count', 'rb'))
    except:
        print(file_kmerT, 'binary file not exist, create new ones')
        npKmer_t, npCount_t = kmerfile2npArrThread(file_kmerfa=file_kmerT, threads=threads, outfile=True,tempFolderInMem=True)
    
    #find target kmers in npKmer_1 npKmer_2
    locs1 = findIdenticalSitesOfQueryInSubject(npKmer_t, npKmer_1)
    locs2 = findIdenticalSitesOfQueryInSubject(npKmer_t, npKmer_2)
    npKmer_1 = npKmer_1[locs1]
    npCount_1 = npCount_1[locs1] / npCount_1.sum()
    npKmer_2 = npKmer_2[locs2]
    npCount_2 = npCount_2[locs2] / npCount_2.sum()
    del npKmer_t, npCount_t
    del locs1, locs2
    file_using = file_kmerT+'.using.binary'
    pickle.dump([npKmer_1, npCount_1, npKmer_2, npCount_2], open(file_using,'wb'))
    del npKmer_1, npCount_1, npKmer_2, npCount_2
    
    
    #split file_fq1, file_fq2 to smaller files
    cmd = 'split {file_fq} -l 400000 -a 5 -d {file_fq}_'.format(file_fq=file_fq1)
    os.system(cmd)
    cmd = 'split {file_fq} -l 400000 -a 5 -d {file_fq}_'.format(file_fq=file_fq2)
    os.system(cmd)
    files_1 = glob.glob('{file_fq}_*'.format(file_fq=file_fq1))
    files_2 = glob.glob('{file_fq}_*'.format(file_fq=file_fq2))
    files_1.sort()
    files_2.sort()
    
    #assign belongs for each pair of seqs
    params = [[fq1, fq2, file_using] for fq1, fq2 in zip(files_1, files_2)]
    pool = Pool(threads)
    results = pool.starmap(splitFq_each, params)
    pool.close()
    print(len(results))
    
    #write files
    os.system('cat {file_fq}_*.0 >{file_fq}.0'.format(file_fq=file_fq1))
    os.system('cat {file_fq}_*.1 >{file_fq}.1'.format(file_fq=file_fq1))
    os.system('cat {file_fq}_*.2 >{file_fq}.2'.format(file_fq=file_fq1))
    os.system('cat {file_fq}_*.0 >{file_fq}.0'.format(file_fq=file_fq2))
    os.system('cat {file_fq}_*.1 >{file_fq}.1'.format(file_fq=file_fq2))
    os.system('cat {file_fq}_*.2 >{file_fq}.2'.format(file_fq=file_fq2))
    
    #remove temp files
    os.system('rm {file_fq}_*'.format(file_fq=file_fq1))
    os.system('rm {file_fq}_*'.format(file_fq=file_fq2))
    print('done')




kmerfile2npArrThread(file_kmer1,threads,True,True)
kmerfile2npArrThread(file_kmer2,threads,True,True)





