import os
from collections import Counter
from Bio import SeqIO
import numpy as np
from multiprocessing import Pool

def npInt8toSeq(baseInt):
    '''
    bases/AAs are coded in a numpy.Int8 array, return a str of bases
    '''
    return ''.join(chr(e) for e in baseInt)

def seq2npInt8(seq):
    '''
    convert bases/AAs to coded in a numpy.Int8 array
    '''
    return np.array([ord(e) for e in seq],dtype=np.int8)

def seqs2npInt8(seqs, threads=16):
    '''
    convert seqs 2 2d array of numpy.int8
    '''
    pool = Pool(threads)
    ls_Int8 = pool.map(seq2npInt8,seqs)
    pool.close()
    return np.array(ls_Int8,dtype=np.int8)

def np2dInt8ToSeqs(np2d,threads=16):
    '''
    convert np2dInt8 to sequences
    '''
    pool = Pool(threads)
    ls_seqs = pool.map(npInt8toSeq,[e for e in np2d])
    pool.close()
    return ls_seqs

def readFasta2np2d(filename,threads=16):
    '''
    filename is a fasta file or a opened file, return a list of sequence names, and a np.ndarray of sequences coding in np.int8
    '''
    ls_seqs = list(SeqIO.parse(filename,'fasta'))
    ls_names = [e.id for e in ls_seqs]
    ls_seqs = [str(e.seq) for e in ls_seqs]
    npInt8seqs = seqs2npInt8(ls_seqs,threads=threads)
    return ls_names,npInt8seqs
def filterSites(np2d, baseKeepMinCount, maxGapCount):
    '''
    np2d is a alignment same as npInt8seqs, alignment coding in np.ndarray
    filter each site based on baseKeepMinCount and maxGapCount
    '''
    keep = []#save sites to keep
    sites_changed = 0
    sites_initial_removed = 0
    sites_further_removed = 0
    sites_only_onebase = 0
    seqNum, seqlen = np2d.shape
    for n in range(seqlen):
        a = np2d[:,n]
        baseCount = Counter(a)
        if baseCount[45]<maxGapCount:
            baseCountKeep = {k:v for k,v in baseCount.items() if v >= baseKeepMinCount}
            if len(baseCount) != len(baseCountKeep):
                sites_changed += 1
            for m in range(seqNum):#change bases less than baseKeepMinCount to gap
                if np2d[m,n] not in baseCountKeep:
                    np2d[m,n] = 45
            a = np2d[:,n]
            baseCount = Counter(a)
            if baseCount[45]<maxGapCount:
                if 45 in baseCount:
                    if len(baseCount) > 2:
                        keep.append(n)
                    else:
                        sites_only_onebase += 1
                else:
                    if len(baseCount) > 1:
                        keep.append(n)
                    else:
                        sites_only_onebase += 1
            else:
                sites_further_removed += 1
        else:
            sites_initial_removed += 1
    print(len(keep), 'sites keep from totally', seqlen,'sites. initial remove', sites_initial_removed,'sites, then ', sites_changed, 'sites were changed.', sites_further_removed,'were further removed based on gap.', sites_only_onebase, 'sites removed as there is only one kind of bases')
    return np2d[:,keep]


#threads = 48
#baseKeepMinCount = 10 #for each position, if the existence of a base is less than baseKeepMinCount, it will be considered as a gap, change that value to gap, which is 45
#minNonGapRatio = 0.8 #after change bases with low frequency to gap, only keep sites with bases ratio grater than minNonGapRatio
#filename = '/work/archive/biophysics/Nick_lab/shared/Xiaolong/20181010Tree/20181010junoniaSeq_z'
#outfile = '/work/archive/biophysics/Nick_lab/shared/Xiaolong/20181010Tree/20181010junoniaSeq_z.filter'

def fastaAlignmentFilter(filename, outfile=None, threads = 16, baseKeepMinCount = 10, minNonGapRatio = 0.8):
    '''
    filename is a input of fasta alignment, output a file filter baseKeepMinCount and minNonGapRatio
    if outfile is None, outfile = filename+'.filter'
    '''
    if outfile is None:
        outfile = filename+'.filter'
    ls_names, npInt8seqs = readFasta2np2d(filename, threads=threads)
    # calculate max allowed gaps per site
    maxGapCount = (1-minNonGapRatio) * npInt8seqs.shape[0]
    # split npInt8seqs to threads parts to use multiprocessing
    seqlen = npInt8seqs.shape[1]
    step = seqlen // threads + 1
    ls_npInt8seqs = [npInt8seqs[:,i:i+step] for i in range(0,seqlen,step)]
    # filter sites
    pool = Pool(threads)
    ls_npInt8seqsKeep = pool.starmap(filterSites,[[e, baseKeepMinCount, maxGapCount] for e in ls_npInt8seqs])
    pool.close()
    npInt8seqsKeep = np.concatenate(ls_npInt8seqsKeep,axis=1)
    ls_seqsKeep = np2dInt8ToSeqs(npInt8seqsKeep, threads=threads)
    
    fout = open(outfile,'w')
    for seq_id, seq in zip(ls_names,ls_seqsKeep):
        fout.write('>'+seq_id+'\n'+seq+'\n')
    fout.close()
    print('done',filename)
    return None

description = '''filename is a input of fasta alignment, output a file filter baseKeepMinCount and minNonGapRatio. Also, only keep sites with information, which means that there should be more than 1 ind of bases per site
    if outfile is None, outfile = filename+'.filter'
default
#threads = 16
#baseKeepMinCount = 10 #for each position, if the existence of a base is less than baseKeepMinCount, it will be considered as a gap, change that value to gap, which is 45
#minNonGapRatio = 0.8 #after change bases with low frequency to gap, only keep sites with bases ratio grater than minNonGapRatio
'''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file of fasta file', required=True)
    parser.add_argument('-o','--outfile', help = 'outfile, where to store the output file. default=None', default=None)
    parser.add_argument('-t','--threads', help = 'threads, number of CPUs to use. default=16', default=16,type=int)
    parser.add_argument('-b','--baseKeepMinCount', help = 'baseKeepMinCount, for each site, keep the base if the base cout is greater than baseKeepMinCount. default=10', default=10,type=int)
    parser.add_argument('-g','--minNonGapRatio', help = 'minNonGapRatio, keep one site only if the number of non-gap bases is greater than this ratio. default = 0.8', default=0.8,type=float)
    f = parser.parse_args()
    fastaAlignmentFilter(filename=f.input, outfile=f.outfile, threads = f.threads, baseKeepMinCount = f.baseKeepMinCount, minNonGapRatio = f.minNonGapRatio)






