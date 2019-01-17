# -*- coding: utf-8 -*-
"""
Created on Sat Jan 12 23:03:54 2019

@author: ATPs
"""
from Bio import SeqIO
import os
from collections import Counter
import pandas as pd
from multiprocessing import Pool

def onlyATCG(seq):
    '''
    all capital letters in seq
    return True if seq only contain ATCG
    '''
    for s in seq:
        if s not in 'ATCG':
            #print(s,'not in ATCG for', seq)
            return False
    return True

def getKmers(filename,dc_kmer2id):
    '''
    a helper function
    '''
    ls_kmers = [e for e in SeqIO.parse(filename,'fasta') if str(e.seq) in dc_kmer2id]
    ls_kmers = ['>'+e.id+'\n'+str(e.seq)+'\n' for e in ls_kmers]
    return ''.join(ls_kmers)


#generate file 
#file_barcodes = '/work/biophysics/s185491/20181226Calephelis/20190110BarcodeChecking/allBarcodes/intermediate_barcodes_20181231.fa'
#file_inputseqs = '/work/biophysics/s185491/20181226Calephelis/20190110BarcodeChecking/fa/3505'
#JELLYFISH = '/home2/s185491/p/anaconda3/anaconda520/envs/bio/bin/jellyfish'
#kmerlen = 19
#threads = 32
#target_sample = '3505'
#outprefix = '/work/biophysics/s185491/20181226Calephelis/20190110BarcodeChecking/fa/3505'

def check_barcode(file_barcodes, file_inputseqs, JELLYFISH, kmerlen,threads, target_sample, outprefix):
    '''
    find the best barcodes for input fasta/fastq file(s). and find possible contaminant barcodes
    file_barcodes: fasta store the barcode sequences
    file_inputseqs: file store the input sequencing seqs. files can be fasta or fastq format. separate by ';' if there are more than one files
    JELLYFISH is where the program jelleyfish stored. The default works on BioHPC
    kmerlen: kmer length to count. default 21
    threads: number of CPUs to use. default 32
    target_sample: the target sample of input sequencing data. default None, program will find the best target barcode by counting kmers. It can help find contaminant target. kmers map to target_samples will be removed when try to find contaminant target barcode
    outprefix: prefix to store the outputs. default is the same as the first elements of file_inputseqs
    '''
    file_inputseqs = file_inputseqs.split(';')
    
    if outprefix is None:
        outprefix = file_inputseqs[0]
    
    file_inputseqs = ' '.join(file_inputseqs)
    
    ls_barcodes = list(SeqIO.parse(file_barcodes,'fasta'))
    ls_ids = [e.id for e in ls_barcodes]
    sample_ids = [e.split('_')[0] for e in ls_ids]
    ls_barseqs = [str(e.seq).upper() for e in ls_barcodes] #keep seqs and change all to capital letters
    print('total barcode seqs',len(ls_barseqs))
    dc_kmer2id = {}
    for n,seq in enumerate(ls_barseqs):
        for i in range(len(seq)+1-kmerlen):
            kmernum = seq[i:i+kmerlen]
            if kmernum not in dc_kmer2id:
                dc_kmer2id[kmernum] = set()
            dc_kmer2id[kmernum].add(n)
    
    # change values to list. remove kmers with none ATCG letters
    
    
    dc_kmer2id = {k:list(v) for k,v in dc_kmer2id.items() if onlyATCG(k)}
    print('total kmers with length',kmerlen,'is', len(dc_kmer2id))
    
    cmd = '{JELLYFISH} count -m {kmerlen} -s 100M -L 5 -t {threads} -C {file_inputseqs} -o {file_inputseqs}.mer_counts.jf && {JELLYFISH} dump {file_inputseqs}.mer_counts.jf >{file_inputseqs}.mer_counts.fa && rm {file_inputseqs}.mer_counts.jf'.format(JELLYFISH=JELLYFISH, kmerlen=kmerlen,threads=threads, file_inputseqs=file_inputseqs)
    os.system(cmd)
    
    #split the file to small files, each with 10,000,000 lines
    file_kmerfa = '{file_inputseqs}.mer_counts.fa'.format(file_inputseqs=file_inputseqs)
    cmd = 'split {file_kmerfa} -l 10000000 -a 5 -d {file_kmerfa}_'.format(file_kmerfa=file_kmerfa)
    os.system(cmd)
    
    import glob
    files = glob.glob('{file_kmerfa}_*'.format(file_kmerfa=file_kmerfa))
    
    #extract fa that exist in dc_kmer2id
    
    pool = Pool(threads)
    results = pool.starmap(getKmers, [[f, dc_kmer2id] for f in files])
    pool.close()
    os.system('rm {file_kmerfa}*'.format(file_kmerfa=file_kmerfa))
    
    file_kmerBar = '{file_inputseqs}.kmer.bar.fa'.format(file_inputseqs=file_inputseqs)
    open(file_kmerBar,'w').write(''.join(results))
    
    ls_kmer_bait = list(SeqIO.parse(file_kmerBar,'fasta'))
    ls_kmer_bait.sort(key=lambda e:int(e.id), reverse=True)
    ls_kmer_bait_counts = [int(e.id) for e in ls_kmer_bait]
    if len(ls_kmer_bait) > 5:
        print('most abundant kmers are')
        for i in range(5):
            print(str(ls_kmer_bait[i].seq), ls_kmer_bait[i].id)
        for i in range(4):
            if ls_kmer_bait_counts[i] / ls_kmer_bait_counts[4] > 3:
                print(i, ls_kmer_bait[i].seq, 'removed. count more than 3 times of the 4th one (starting from 0)')
            else:
                break
        ls_kmer_bait = ls_kmer_bait[i:]
        ls_kmer_bait_counts = ls_kmer_bait_counts[i:]
        
        
    total_unique_kmer = len(ls_kmer_bait)
    print('total unique kmers', total_unique_kmer)
    total_kmer = sum(int(e.id) for e in ls_kmer_bait)
    print('total kmers matched to bait', total_kmer)
    
    #count kmer for each input
    barcodeCounter = Counter()
    barcodeCounter_unique = Counter()
    for kmer in ls_kmer_bait:
        seq = str(kmer.seq)
        count = int(kmer.id)
        kmerCounter = Counter(dc_kmer2id[seq])
        barcodeCounter_unique = barcodeCounter_unique+kmerCounter
        for k in kmerCounter:
            kmerCounter[k] = kmerCounter[k] * count
        barcodeCounter = barcodeCounter + kmerCounter
    
    #get best id
    df_result = pd.DataFrame(barcodeCounter.most_common())
    df_result.columns = ['ID_in_barcodefile','kmer_identified']
    df_result = df_result.set_index('ID_in_barcodefile')
    df_result['sample_id'] = [ls_ids[e] for e in df_result.index]
    df_result['total_kmer'] = total_kmer
    df_result['total_unique_kmer'] = total_unique_kmer
    df_result['unique_kmer_identified'] = [barcodeCounter_unique[e] for e in df_result.index]
    df_result['kmer_ratio'] = df_result['kmer_identified'] / df_result['total_kmer']
    df_result['unique_kmer_ratio'] = df_result['unique_kmer_identified'] / df_result['total_unique_kmer']
    df_result = df_result[['sample_id', 'kmer_identified', 'total_kmer', 'kmer_ratio', 'unique_kmer_identified', 'total_unique_kmer', 'unique_kmer_ratio']]
    
    #count kmer that do not belongs to the target_sample
    if target_sample not in sample_ids:
        print(target_sample, 'not in the input barcode. set to None')
        target_sample = None
    if target_sample is None:
        print('target sample is not provided, use the ones with highest kmer_identified, which are')
        target_sample = df_result[df_result['kmer_identified'] == df_result['kmer_identified'].max()]
        for s in target_sample['sample_id']:
            print(s)
    else:
        target_sample = df_result[df_result['sample_id'].apply(lambda x:x.split('_')[0]==target_sample)]
    
    target_id = list(target_sample.index)
    kmer_to_exclude = set()
    for k,v in dc_kmer2id.items():
        if len(v) + len(target_id) != len(set(v+target_id)):
            kmer_to_exclude.add(k)
    ls_kmer_bait = [e for e in ls_kmer_bait if str(e.seq) not in kmer_to_exclude]
    barcodeCounter = Counter()
    barcodeCounter_unique = Counter()
    for kmer in ls_kmer_bait:
        seq = str(kmer.seq)
        count = int(kmer.id)
        kmerCounter = Counter(dc_kmer2id[seq])
        barcodeCounter_unique = barcodeCounter_unique+kmerCounter
        for k in kmerCounter:
            kmerCounter[k] = kmerCounter[k] * count
        barcodeCounter = barcodeCounter + kmerCounter
    df_result['nonTargetKmerCount'] = [barcodeCounter[e] for e in df_result.index]
    df_result['nonTarget_uniqueKmerCount'] = [barcodeCounter_unique[e] for e in df_result.index]
    
    # most possible contamination
    df_result2 = df_result.sort_values(by='nonTargetKmerCount',ascending=False)
    df_result2['nonTargetKmerCount_ratio'] = df_result2['nonTargetKmerCount'] / df_result2['total_kmer']
    
    #write the result
    os.rename(file_kmerBar, outprefix+'.kmer.bar.fa')
    df_result.head(100).to_csv(outprefix+'.kmercount.csv',sep='\t')
    file_reports = outprefix+'.reports'
    fout = open(file_reports,'w')
    reports = '''input file: {file_inputseqs}
    kmerlen used in analysis: {kmerlen}
    total barcode seqs: {barcodeseqs}
    
    best match in barcodes:
    {bestmatch}
    
    target match in barcodes:
    {targetmatch}
    
    most possible contamination:
    {contamination}
    '''.format(file_inputseqs=file_inputseqs, kmerlen=kmerlen, barcodeseqs=len(ls_barseqs),bestmatch = df_result[df_result['kmer_identified'] == df_result['kmer_identified'].max()].to_string(), contamination = df_result2[df_result2['nonTargetKmerCount'] == df_result2['nonTargetKmerCount'].max()].to_string(), targetmatch = target_sample.to_string())
    fout.write(reports)
    fout.close()
    print('done')

description = '''
find the best barcodes for input fasta/fastq file(s). and find possible contaminant barcodes
    file_barcodes: fasta store the barcode sequences
    file_inputseqs: file store the input sequencing seqs. files can be fasta or fastq format. separate by ';' if there are more than one files
    JELLYFISH is where the program jelleyfish stored. The default works on BioHPC
    kmerlen: kmer length to count. default 21
    threads: number of CPUs to use. default 32
    target_sample: the target sample of input sequencing data. default None, program will find the best target barcode by counting kmers. It can help find contaminant target. kmers map to target_samples will be removed when try to find contaminant target barcode. it should be like '3505', the prefix. the program will split('_') and use the first element
    outprefix: prefix to store the outputs. default is the same as the first elements of file_inputseqs
'''
if __name__ == '__main__':
    print(description)
    import argparse
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-b','--file_barcodes', help = 'where the barcodes were stored', required=True)
    parser.add_argument('-i','--file_inputseqs', help = '''file_inputseqs is file store the input sequencing seqs. files can be fasta or fastq format. separate by ';' if there are more than one files''')
    parser.add_argument('-J','--JELLYFISH', help = ''' where the program jelleyfish stored. The default works on BioHPC. default = "/home2/s185491/p/anaconda3/anaconda520/envs/bio/bin/jellyfish". Note, when running jellyfish, only kmers with counts more than 5 will be reported. Top 4 abundant kmers matching the barcodes, if they are more than 3 times of the 4th kmer, they would be excluded from the analysis, because some kmers a just too universal in the sequenced reads ''', default = '/home2/s185491/p/anaconda3/anaconda520/envs/bio/bin/jellyfish')
    parser.add_argument('-k','--kmerlen', help = '''kmer length to count. default 21''', default=21, type=int)
    parser.add_argument('-t','--threads', help = '''number of CPUs to use. default 32''', default=32, type=int)
    parser.add_argument('-T','--target_sample', help = '''the target sample of input sequencing data. default None, program will find the best target barcode by counting kmers. It can help find contaminant target. kmers map to target_samples will be removed when try to find contaminant target barcode,it should be like '3505', the prefix. the program will split('_') and use the first element''', default=None)
    parser.add_argument('-o','--outprefix', help = '''three output file will be generated. the prefix for those files. default is the same as the first elements of file_inputseqs''', default = None)
    f = parser.parse_args()
    check_barcode(file_barcodes=f.file_barcodes, file_inputseqs=f.file_inputseqs, JELLYFISH=f.JELLYFISH, kmerlen=f.kmerlen,threads=f.threads, target_sample=f.target_sample, outprefix=f.outprefix)