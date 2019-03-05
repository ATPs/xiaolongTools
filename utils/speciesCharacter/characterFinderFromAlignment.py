# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:59:18 2019

@author: ATPs
"""

#file_alignment = '/work/archive/biophysics/Nick_lab/jzhang/project/jzhang/Ortholexis/Ortholexis_nuclear_20181216.fasta'
#threads = 32
#file_sampleInfo = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/20190122JingGroup4_sample_info'
#file_loc = '/work/archive/biophysics/Nick_lab/jzhang/project/jzhang/Ortholexis/Ortholexis_dna_range_v1'
#output_prefix = '/work/biophysics/s185491/20181226Calephelis/20190116geneIsland/JingExample/20190122JingGroup4_5'
#file_tree = '/archive/biophysics/Nick_lab/jzhang/project/jzhang/tmp/Ortholexis_nuclear_20181216.newick.rescaleAtt1.5.renamed'
#file_annotation = '/archive/biophysics/Nick_lab/jzhang/project/jzhang/pyrr/aly_table'
#nonGapRatio_target = 1
#nonGapRatio_close=0.9
#dominantRatio_close = 0.9 #ratio calculation with gaps excluded

import os
from collections import Counter
from Bio import SeqIO
import numpy as np
from multiprocessing import Pool
import pandas as pd
import ete3


def npInt8toSeq(baseInt):
    '''
    bases are coded in a numpy.Int8 array, return a str of bases
    '''
    return ''.join(chr(e) for e in baseInt)

def seq2npInt8(seq, N2gap=True):
    '''
    convert bases to coded in a numpy.Int8 array
    N2gap is True, convert "N" to "-"
    '''
    seq = seq.upper()
    if N2gap:
        seq = seq.replace('N','-')
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

def cleanSampleId(sample_id):
    '''
    if sample_id is long description like "18082A08_Ortholexis_hollandi_M_old"
    return "18082A08"
    '''
    return sample_id.split("_")[0]

def readSampleInfo(file_sampleInfo):
    '''
    return a dataframe with sample_id, sampe_description, phylogroup (OCTE) and group. with sample_id as index
    each sample will have a group like T1, T2, or T, O, E
    O: outgroup
    C: close group: typically those groups not the target group, not outgroup
    E: samples to be excluded from the analysis. But will be included in the output alignment file
    T: target samples to study
    '''
    df_sample_info = pd.read_csv(file_sampleInfo,sep='\t', dtype=str)
    df_sample_info = df_sample_info.set_index('sample_id')
    #replace space with '_' for sample_description
    df_sample_info['sample_description'] = df_sample_info['sample_description'].apply(lambda x:x.replace(' ', '_'))
    df_sample_info['phylogroup'] = df_sample_info['group'].apply(lambda x:x[0])
    return df_sample_info
    

def readFastaAlignmentFile2df(file_alignment,threads):
    '''
    file_alignment is a fasta alignment file.
    return a dataframe, with columns of sample_id, and index of positions starting from 0. base stored in np.int8
    '''
    ls_names, npInt8seqs = readFasta2np2d(filename=file_alignment, threads=threads)
    ls_names = [cleanSampleId(e) for e in ls_names]
    df_seqs = pd.DataFrame(npInt8seqs)
    df_seqs = df_seqs.T
    df_seqs.columns = ls_names
    return df_seqs

def combineNpInt8Seq(seqs):
    seqlen = len(seqs[0])
    if len(seqs) == 1:
        return seqs[0]
    results = np.ones(seqlen, dtype=np.int8) * 45
    for n, bases in enumerate(zip(*seqs)):
        bases_nongap = [e for e in bases if e not in [78,45]]
        if len(bases_nongap) != 0:
            if len(set(bases_nongap)) == 1:
                results[n] = bases_nongap[0]
            else:
                results[n] = 78
    return results

def combineNpInt8SeqThreads(seqs, threads):
    '''
    seqs is a list of npInt8 array, with int8 represent ATCGN-
    return a combined seq in npInt8
    A 65
    T 84
    C 67
    G 71
    N 78
    - 45
    if there is only one kind of base, return that base. If only with N and -, return -. else, return N
    '''
    if len(seqs) == 1:
        return seqs[0]
    ls_seqs = [np.array_split(e, threads) for e in seqs]
    ls_seqs = list(zip(*ls_seqs)) #split seqs to threads fragments
    pool = Pool(threads)
    ls_results = pool.map(combineNpInt8Seq, ls_seqs)
    pool.close()
    return np.concatenate(ls_results)

def countBasesNpInt8Seq(seqs, target_bases):
    '''
    seqs is a list of npInt8 sequences, return a npInt16 array with equal size
    target_bases is a list bases coded in int or a string of "ATCGN-"
    '''
    if isinstance(target_bases, str):
        target_bases = [ord(e) for e in target_bases]
    seqlen = len(seqs[0])
    results = np.zeros(seqlen,dtype=np.int16)
    for n, bases in enumerate(zip(*seqs)):
        results[n] = len([e for e in bases if e in target_bases])
    return results

def countBasesNpInt8SeqThreads(seqs, target_bases, threads):
    '''
    seqs is a list of npInt8 sequences, return a npInt16 array with equal size
    use multiple thread to accelarate
    target_bases is a list bases coded in int or a string of "ATCGN-"
    '''
    ls_seqs = [np.array_split(e, threads) for e in seqs]
    ls_seqs = list(zip(*ls_seqs)) #split seqs to threads fragments
    pool = Pool(threads)
    ls_results = pool.starmap(countBasesNpInt8Seq, [[e, target_bases] for e in ls_seqs])
    pool.close()
    return np.concatenate(ls_results)

def countDominantBasesNpInt8Seq(seqs):
    '''
    seqs is a list of npInt8 sequences, return a npInt16 array with equal size
    and a npInt8 array with dominant bases. dominant bases cannot be '-N'. If with only '-N', use '-'
    '''
    seqlen = len(seqs[0])
    results_count = np.zeros(seqlen,dtype=np.int16)
    results_base = np.zeros(seqlen,dtype=np.int8)
    for n, bases in enumerate(zip(*seqs)):
        bases_nogap = [e for e in bases if e not in [78,45]]
        baseCounter = Counter(bases_nogap).most_common()
        if len(baseCounter) == 0:
            results_count[n] = 0
            results_base[n] = 45
        else:
            results_count[n] = baseCounter[0][1]
            results_base[n] = baseCounter[0][0]
    return results_base, results_count

def countDominantNpInt8SeqThreads(seqs, threads):
    '''
    seqs is a list of npInt8 sequences, return a npInt16 array with equal size
    use multiple thread to accelarate
    target_bases is a list bases coded in int or a string of "ATCGN-"
    '''
    ls_seqs = [np.array_split(e, threads) for e in seqs]
    ls_seqs = list(zip(*ls_seqs)) #split seqs to threads fragments
    pool = Pool(threads)
    ls_results = pool.map(countDominantBasesNpInt8Seq, ls_seqs)
    pool.close()
    ls_results_base = [e[0] for e in ls_results]
    ls_results_count = [e[1] for e in ls_results]
    return np.concatenate(ls_results_base), np.concatenate(ls_results_count)

def getWorkingSeqs(df_seqs, df_sample_info, phylogroups, threads, combineGroup=True):
    '''
    return a list of seqs in npInt8 format that belongs to phylogroups.
    phylogroups is a string with "OTCE"
    if combineGroup is True, combine samples based on group in df_sample_info. 
    Note: single letter means single sample in that group. only combine group like T1, T0
    '''
    df_sample_keep = df_sample_info[df_sample_info['phylogroup'].apply(lambda x: x in phylogroups)].copy()
    if not combineGroup:
        sample_ids_keep = df_sample_keep.index
        ls_seqs = [np.array(df_seqs[e]) for e in sample_ids_keep]
    else:
        sample_ids_keep = df_sample_keep[df_sample_keep['group'].apply(lambda x:len(x)==1)].index
        ls_seqs = [np.array(df_seqs[e]) for e in sample_ids_keep]
        df_sample_toCombine = df_sample_keep[df_sample_keep['group'].apply(lambda x:len(x)!=1)]
        for group_toCombine in df_sample_toCombine['group'].unique():
            samples_toCombine = df_sample_toCombine[df_sample_toCombine['group']==group_toCombine].index
            seqs_toCombine = [np.array(df_seqs[e]) for e in samples_toCombine]
            ls_seqs.append(combineNpInt8SeqThreads(seqs_toCombine, threads))
    return ls_seqs


def readFragLoc(file_loc, exonSpliter='.e'):
    '''
    file_loc is a file of exon location of genes
    each line looks like
        #lac1000.1.e2	0	69
    exonSpliter is used to rsplit the exon_id to get gene_id and exon_N
    if exonSpliter is None, gene_id = exon_id
    return a dataframe
    '''
    df = pd.read_csv(file_loc,sep='\t',header=None)
    df.columns = ['exon_id','start','end']
    if exonSpliter is None or exonSpliter == 'None':
        df['gene_id'] = df['exon_id']
        df['exon_N'] = 1
    else:
        df['gene_id'] = df['exon_id'].apply(lambda x:x.rsplit(exonSpliter,1)[0])
        df['exon_N'] = df['exon_id'].apply(lambda x:int(x.rsplit(exonSpliter,1)[1]))
    df = df.sort_values(by = ['gene_id','exon_N'])
    df['exon_len'] = df['end'] - df['start']
    df['exon_start_gene'] = df.groupby('gene_id')['exon_len'].cumsum() - df['exon_len']
    df['gene_len'] = df['exon_len'].groupby(df['gene_id']).transform('sum')
    return df

def getExonGeneInfoForSite(sites, df_loc):
    '''
    df_loc looks like:
              exon_id     start       end  gene_id  exon_N  exon_len  exon_start_gene  gene_len
    11862   aly1.1.e1   1885263   1885389   aly1.1       1       126                0       729
    29574   aly1.1.e2   5657499   5658102   aly1.1       2       603              126       729
    78034   aly1.2.e1  16756914  16759644   aly1.2       1      2730                0      8280
    72098   aly1.2.e2  15367371  15370080   aly1.2       2      2709             2730      8280
    79838   aly1.2.e3  17197977  17200818   aly1.2       3      2841             5439      8280
    return exon_id, gene_id, location in exon of site, location of gene
    eg, sites = [5657500], return [['aly1.1.e2', 'aly1.1', 1, 127]]
    '''
    df = df_loc.copy()
    df = df.sort_values(by = 'start')
    start_array = np.array(df['start'])
    end_array = np.array(df['end'])
    index_basedonStart = np.searchsorted(start_array, sites, side='right')-1
    index_basedonEnd = np.searchsorted(end_array, sites, side='right')
    if not all(index_basedonStart == index_basedonEnd):
        print('something wrong with the file_loc, exon location overlap or other problem')
        return None
    df_keep = df.iloc[index_basedonEnd].copy()
    df_keep['sites'] = sites
    df_keep['loc_exon'] = df_keep['sites'] - df_keep['start']
    df_keep['loc_gene'] = df_keep['loc_exon'] + df_keep['exon_start_gene']
    df_keep = df_keep.set_index('sites')
    return df_keep


def generateSNPdescription(row, loc_column='loc_exon'):
    '''
    generate description like G3092A, G is dominantNt_close, A is dominantNt_target, 3092 is the location
    if loc_column is None, return GtoA
    '''
    nt_other = row['dominantNt_close']
    nt_target = row['dominantNt_target']
    if loc_column is None:
        loc = 'to'
    else:
        loc = str(row[loc_column])
    return nt_other+loc+nt_target

def getSampleNamesFromTreeFile(file_tree):
    '''
    return a list of sample_ids, based the tree file
    '''
    tree = ete3.Tree(file_tree)
    names = [node.name for node in tree]
    return [name.strip("'").strip('"').split('_')[0] for name in names]


def writeAlignment(df_seqs, df_simple, df_sample_info, names_order, outfilename, fraglen = 21):
    '''
    df_seqs is sequences in dataframe
    df_simple is description about SNPs
    names_order is sample_ids to use
    df_sample_info contails the long description for each sample
    outfilename is where to store the output
    fraglen is the length of fragment. In most time, SNP site is in the middle of the fragments, unless the site too close to the end
    '''
    fout = open(outfilename,'w')
    N = max(df_sample_info['sample_description'].apply(len))+1
    for n, row in df_simple.iterrows():
        try:
            txt = '###alignment for {gene_id} {exon_id} exon: {SNPloc_exon} (gapClose: {gapCount_close}; SNPs: {gene_SNPcount}) ###\n'.format(gene_id=row['gene_id'], exon_id=row['exon_id'], SNPloc_exon=row['SNPloc_exon'], gapCount_close=row['gapCount_close'], gene_SNPcount=row['gene_SNPcount'])
        except:
            txt = '###alignment for {SNPloc_alignment} (gapClose: {gapCount_close}) ###\n'.format(SNPloc_alignment=row['SNPloc_alignment'], gapCount_close=row['gapCount_close'])
        fout.write(txt)
        start = max(n - (fraglen//2),0)
        end = start + fraglen
        df_seqs_frag = df_seqs.loc[list(range(start,end))].copy()
        dc_seq_frag = df_seqs_frag.apply(lambda x:''.join(chr(e) for e in x), axis='index').to_dict()
        df_frag = df_sample_info.loc[names_order,['sample_description']].copy()
        df_frag['seq'] = [dc_seq_frag[e] for e in df_frag.index]
        for nn, rrow in df_frag.iterrows():
            fout.write('{description:<{N}} {seq}\n'.format(N=N, description=rrow['sample_description'], seq=rrow['seq']))
        fout.write('\n')
    fout.close()
    return None


def readFileAnnotation(file_annotation):
    '''
    file_annotation is a annotation file for genes. 
    the first column is gene name and the second column is annotation for genes
    '''
    dc_gene2anno = {}
    for line in open(file_annotation):
        items = line.strip().split('\t')
        dc_gene2anno[items[0]] = items[1]
    return dc_gene2anno



def characterFinder(file_alignment, file_sampleInfo, file_loc=None, file_annotation=None, file_tree=None,nonGapRatio_target = 1, nonGapRatio_close=0.9, dominantRatio_target = 1, dominantRatio_close = 0.9, threads=32, output_prefix=None,exonSpliter='.e'):
    '''
    '''
    # read in sample info
    df_sample_info = readSampleInfo(file_sampleInfo)
    samples_exclude = list(df_sample_info[df_sample_info['phylogroup'] == 'E'].index)
    samples_target = list(df_sample_info[df_sample_info['phylogroup'] == 'T'].index)
    samples_close = list(df_sample_info[df_sample_info['phylogroup'] == 'C'].index)
    samples_outgroup = list(df_sample_info[df_sample_info['phylogroup'] == 'O'].index)
    samples_use = samples_target + samples_close +samples_outgroup
    print('finished reading in file_sampleInfo. Based on this file, there are', df_sample_info.shape[0], 'samples')
    print('target samples', samples_target)
    print('outgroup samples', samples_outgroup)
    print('close group samples', samples_close)
    print('samples to exclude', samples_exclude)
    print(len(samples_use), 'samples will be used for analysis')
    
    # read in fasta alignment file
    df_seqs = readFastaAlignmentFile2df(file_alignment,threads)
    print("finish reading in alignment.")
    print("sequence length", df_seqs.shape[0])
    print("number of seqs in the alignment file", df_seqs.shape[1])
    
    # remove seqs that have a phylogroup of 'E'
    #df_seqs = df_seqs.drop(labels = samples_exclude, axis='columns')
    print(len(samples_exclude),'''samples labeled with "E" and not included from analysis''')
    
    #count gaps for target_group
    ls_seqs_target = getWorkingSeqs(df_seqs=df_seqs, df_sample_info=df_sample_info, phylogroups = 'T', threads=threads, combineGroup=True)
    seqCount_target = len(ls_seqs_target)
    print(len(samples_target),'samples in target group, after combine samples belong to same cluster, the remaining seqs are', seqCount_target)
    npGaps_target = countBasesNpInt8SeqThreads(seqs=ls_seqs_target, target_bases='-', threads=threads)
    df_seqs_filter = df_seqs[npGaps_target <= seqCount_target * (1 - nonGapRatio_target)].copy()
    print('keep sites that the target samples have nonGapRatio >=', nonGapRatio_target)
    print(df_seqs_filter.shape[0],'sites left. {percent:.2%} of total'.format(percent = df_seqs_filter.shape[0] / df_seqs.shape[0]))
    
    #count gaps for close_group
    ls_seqs_close = getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'C', threads=threads, combineGroup=True)
    seqCount_close = len(ls_seqs_close)
    print(len(samples_close),'samples in target group, after combine samples belong to same cluster, the remaining seqs are', seqCount_close)
    npGaps_close = countBasesNpInt8SeqThreads(seqs=ls_seqs_close, target_bases='-', threads=threads)
    df_seqs_filter = df_seqs_filter[npGaps_close <= seqCount_close * (1 - nonGapRatio_close)].copy()
    print('keep sites that the close samples have nonGapRatio >=', nonGapRatio_close)
    print(df_seqs_filter.shape[0],'sites left. {percent:.2%} of total'.format(percent = df_seqs_filter.shape[0] / df_seqs.shape[0]))
    
    #count dominant bases for target_group and close_group
    ls_seqs_target = getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'T', threads=threads, combineGroup=True)
    ls_seqs_close = getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'C', threads=threads, combineGroup=True)
    npDominantBase_target, npDominantCount_target = countDominantNpInt8SeqThreads(seqs=ls_seqs_target, threads=threads)
    npDominantBase_close, npDominantCount_close = countDominantNpInt8SeqThreads(seqs=ls_seqs_close, threads=threads)
    df_seqs_filter = df_seqs_filter[npDominantBase_target != npDominantBase_close]
    print('dominant bases are not the same in target group and closegroup')
    print(df_seqs_filter.shape[0],'sites left. {percent:.2%} of total'.format(percent = df_seqs_filter.shape[0] / df_seqs.shape[0]))
    
    ls_seqs_target = getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'T', threads=threads, combineGroup=True)
    ls_seqs_close = getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'C', threads=threads, combineGroup=True)
    npDominantBase_target, npDominantCount_target = countDominantNpInt8SeqThreads(seqs=ls_seqs_target, threads=threads)
    npDominantBase_close, npDominantCount_close = countDominantNpInt8SeqThreads(seqs=ls_seqs_close, threads=threads)
    df_result = df_seqs_filter.copy()
    df_result['seqCount_target'] = seqCount_target
    df_result['seqCount_close'] = seqCount_close
    df_result = df_result.drop(labels = samples_use, axis='columns')
    df_result['gapCount_target'] = countBasesNpInt8SeqThreads(seqs=ls_seqs_target, target_bases='-', threads=threads)
    df_result['gapCount_close'] = countBasesNpInt8SeqThreads(seqs=ls_seqs_close, target_bases='-', threads=threads)
    df_result['dominantBase_target'] = npDominantBase_target
    df_result['dominantCount_target'] = npDominantCount_target
    df_result['dominantBase_close'] = npDominantBase_close
    df_result['dominantCount_close'] = npDominantCount_close
    
    df_result['dominantRatio_target'] = df_result['dominantCount_target'] / (df_result['seqCount_target'] - df_result['gapCount_target'])
    df_result = df_result[df_result['dominantRatio_target'] >= dominantRatio_target]
    print('keep those sites that domiant base greater than {dominantRatio_target} in target samples (exclude gap sites), {n} sites left'.format(dominantRatio_target=dominantRatio_target, n = df_result.shape[0]))
    
    df_result['dominantRatio_close'] = df_result['dominantCount_close'] / (df_result['seqCount_close'] - df_result['gapCount_close'])
    df_result =df_result[df_result['dominantRatio_close'] >= dominantRatio_close]
    print('keep those sites that domiant base ratio (exclude gap sites) greater than', dominantRatio_close, 'in close group samples')
    print('finally {n} sites left'.format(n = df_result.shape[0]))
    
    if df_result.shape[0] == 0:
        print('nothing found, no good. change the parameters!')
        return None
    
    df_result['seqCount_target_ori'] = len(samples_target)
    df_result['seqCount_close_ori'] = len(samples_close)
    df_seqs_filter = df_seqs.loc[df_result.index]
    df_result['gapCount_target_ori'] = countBasesNpInt8SeqThreads(seqs=getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'T', threads=threads, combineGroup=False), target_bases='-', threads=threads)
    df_result['gapCount_close_ori'] = countBasesNpInt8SeqThreads(seqs=getWorkingSeqs(df_seqs=df_seqs_filter, df_sample_info=df_sample_info, phylogroups = 'C', threads=threads, combineGroup=False), target_bases='-', threads=threads)
    
    df_result['dominantNt_target'] = df_result['dominantBase_target'].apply(lambda x:chr(x))
    df_result['dominantNt_close'] = df_result['dominantBase_close'].apply(lambda x:chr(x))
    df_result['nonConservedNt_close'] = df_result['seqCount_close'] - df_result['gapCount_close'] - df_result['dominantCount_close'] 
    # read in file_loc with location of gene/exon/other fragments
    if file_loc is None:
        print('no file provided with the location of exons/genes')
        df_allInfo = df_result.copy()
        df_allInfo['loc_alignment'] = df_allInfo.index
        df_allInfo['loc_alignment'] = df_allInfo.index
        df_allInfo['SNPloc_alignment'] = df_allInfo.apply(lambda x: generateSNPdescription(x,'loc_alignment'), axis=1)
        df_allInfo = df_allInfo.sort_values(by=['gapCount_target_ori','gapCount_close_ori'],ascending=[True,True])
        df_simple = df_allInfo[['SNPloc_alignment', 'gapCount_target', 'gapCount_close','nonConservedNt_close']].copy()
        
    else:
        df_loc = readFragLoc(file_loc, exonSpliter=exonSpliter)
        print('finished readin exon location')
    
        # add gene_id, exon_id, location in exon
        sites = list(df_result.index)
        df_sites_info = getExonGeneInfoForSite(sites=sites, df_loc=df_loc)
        df_allInfo = pd.merge(left=df_result, right=df_sites_info, left_index=True, right_index=True, how='inner')
        df_allInfo['gene_SNPcount'] = df_allInfo.groupby('gene_id')['gene_id'].transform('count')
        df_allInfo['SNPloc_exon'] = df_allInfo.apply(lambda x: generateSNPdescription(x,'loc_exon'), axis=1)
        df_allInfo['SNPloc_gene'] = df_allInfo.apply(lambda x: generateSNPdescription(x,'loc_gene'), axis=1)
        df_allInfo['loc_alignment'] = df_allInfo.index
        df_allInfo['SNPloc_alignment'] = df_allInfo.apply(lambda x: generateSNPdescription(x,'loc_alignment'), axis=1)
        df_allInfo = df_allInfo.sort_values(by=['gene_SNPcount','gene_id','gapCount_target_ori','gapCount_close_ori'],ascending=[False,True,True,True])
        df_simple = df_allInfo[['gene_id','exon_id','gene_SNPcount','SNPloc_exon', 'gapCount_target', 'gapCount_close','nonConservedNt_close']].copy()
        
    
    if output_prefix is None:
        output_prefix = file_alignment
    output_allInfo = output_prefix + '.allInfo'
    df_allInfo.to_csv(output_allInfo, sep='\t',index=None)
    output_simpleInfo = output_prefix+'.simpleInfo'
    df_simple.to_csv(output_simpleInfo, sep='\t',index=None)
    print('finish writing allInfo and .simpleInfo')
    
    #write alignment with target_group first
    samples_all = samples_target + samples_close +samples_outgroup + samples_exclude
    outfile_alignOrdered = output_prefix + '.alignment.ordered'
    writeAlignment(df_seqs=df_seqs, df_simple=df_simple, df_sample_info=df_sample_info, names_order=samples_all, outfilename=outfile_alignOrdered, fraglen = 21)
    
    #if tree provided, write alignment with tree_order
    if file_tree is not None:
        print('file_tree is provided, output alignment based on tree')
        samples_orderedByTree = getSampleNamesFromTreeFile(file_tree=file_tree)
        outfile_orderedByTree = output_prefix + '.alignment.tree'
        writeAlignment(df_seqs=df_seqs, df_simple=df_simple, df_sample_info=df_sample_info, names_order=samples_orderedByTree, outfilename=outfile_orderedByTree, fraglen = 21)
    
    #if gene annotation file provided, write a file with gene annotation
    if file_annotation is not None and file_loc is not None:
        print('file_annotation provided, write a table with gene annotation')
        dc_annotation = readFileAnnotation(file_annotation=file_annotation)
        from collections import defaultdict
        dc_annotation = defaultdict(None, dc_annotation)
        df_simple['gene_annotation'] = df_simple['gene_id'].apply(lambda x:dc_annotation[x])
        output_simpleAnno = output_prefix+'.simpleAnnotation'
        df_simple.to_csv(output_simpleAnno, sep='\t',index=None)
    
    print('done!')



description = '''
a function to find the conserved letter in sequence alignment. The inputs are:
    
    file_alignment: sequence alignment file
    
    file_sampleInfo: information about samples. It looks like:
        
        sample_id	sample_description	group
        7865	7865_Bibasis mahintha_F_Myanmar_2001-09-29	O
        5270	5270_Burara striata_F_China_Sichuan Prov._2015-08-11	O
        17119G11	17119G11_Hasora __Parata_ chromus_2002-05-30	O
        17118G02	17118G02_Hasora badra badra_F_1989-04-22	O
        
        Note: header lines are required. Headers should be the same, the order of headers can be changed.
        It have three columns, first column is sample_id, second is long description for each sample, the last one is the group of samples.
        
        O: outgroup
        C: close group: typically those groups not the target group, not outgroup
        E: samples to be excluded from the analysis. But will be included in the output alignment file
        T: target samples to study
        each sample will have a group like T1, T2, or T, O, E, C. if more than one samples were noted with T1, C1, samples belong to the same group will be combined
    
    file_loc. default None, do not assign SNP to genes/exons
        file_loc is a file of exon location of genes
        each line looks like
            #lac1000.1.e2	0	69
    exonSpliter is used to rsplit the exon_id to get gene_id and exon_N. default = '.e'. the function uses 
        if exonSpliter is None, gene_id = exon_id
    
    file_annotation: if None,do not output annotation file
    
    file_tree: if None, do not output alignment based on samples in the tree
    
    nonGapRatio_target default 1, 
    
    nonGapRatio_close default 0.9, 
    
    dominantRatio_target default 1,
    
    dominantRatio_close default 0.9, 
    
    threads: default 32, numbers of CPUs to use
    
    output_prefix: default None, the same as file_alignment
'''

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='find the feature letter in sequence alignment')
    print(description)
    parser.add_argument('-a','--file_alignment', help = 'fasta alignment of all samples', required=True)
    parser.add_argument('-i','--file_sampleInfo', help = 'a file of sampl_id and long description for that sample, and group of each sample', required=True)
    parser.add_argument('-l','--file_loc', help = 'a file with the locatio of exons in the fasta alignment. default=None', default=None)
    parser.add_argument('-A','--file_annotation', help = 'a file with protein_id and annotation. typically a result of protein function annotation. default=None', default=None)
    parser.add_argument('-t','--file_tree', help = 'a tree file which will be used to re-order samples, default=None', default=None)
    parser.add_argument('-S','--nonGapRatio_target', help = 'percent of letters in target samples. a float number to filter the gap ratio of InGroup samples. default 1.0', type=float, default=1.0)
    parser.add_argument('-C','--nonGapRatio_close', help = 'percent of letters in closegroup. a float number to filter the gap ratio of close samples. default 0.9', type=float, default=0.9)
    parser.add_argument('-M','--dominantRatio_target', help = 'the ratio of dominant bases in target samples (exclude gaps when calculate ratio). default 1', type=float, default=1)
    parser.add_argument('-D','--dominantRatio_close', help = 'the ratio of dominant bases in close samples (exclude gaps when calculate ratio). default 0.9', type=float, default=0.9)
    parser.add_argument('-o','--output_prefix', help = '''prefix for the result outfiles. default None, which will be file_alignment. default=None''', default=None)
    parser.add_argument('-T','--threads', help = '''number of CPUs to use. default 32''', default=32, type=int)
    parser.add_argument('-E','--exonSpliter', help='''spliter to use when convert exon_id to gene_id. default '.e', which 'aly1.1.e1' gene_id is 'aly1.1'. use rsplit''', default='.e')
    f = parser.parse_args()
    characterFinder(file_alignment=f.file_alignment, file_sampleInfo=f.file_sampleInfo, file_loc=f.file_loc, file_annotation=f.file_annotation, file_tree=f.file_tree, nonGapRatio_target = f.nonGapRatio_target, nonGapRatio_close=f.nonGapRatio_close, dominantRatio_target=f.dominantRatio_target, dominantRatio_close = f.dominantRatio_close, threads=f.threads,  output_prefix=f.output_prefix, exonSpliter=f.exonSpliter)
