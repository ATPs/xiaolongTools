# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 22:14:59 2018

@author: ATPs

map2treemix
"""

import os
from collections import Counter
from Bio import SeqIO
import pandas as pd
import numpy as np

#file_groups = '/home/xcao/w/20180905Junonia_coenia/20181010Trees/20181021treemix/263sample_pop_cdsTree19groups'
sitesUse = False

def _getMaxAllowedGaps(ls_maps, allStrand=True, gapCut=0.8):
    '''
    ls_gaps is a list of mapfile locations. return the maximum allowed gaps
    '''
    ls_mapOpen = [open(e) for e in ls_maps]
    for lines in zip(*ls_mapOpen):
        break
    for f in ls_mapOpen:
        f.close()
    bases = []
    for line in lines:
        es = line.split()
        if allStrand:
            bases += es
        else:
            bases.append(es[0])
    return len(bases) * (1-gapCut)

def _processOneLine(lines, allStrand=True,max_gap=float('+inf'),Ncut=0):
    '''
    
    '''

def maps2treemixInput(maps, allStrand=True, siteKeep=None, Ncut = 4, gapCut=0.8, output=None, threads=16):
    '''
    maps is a file of mapfile locations and their species name, each line looks like
        'filename  species'
    maps can also be a list, each element with two parts: filename and species
    siteKeep is a text file which determines sites in the mapfile to keep.
        default, siteKeep=None, all sites will be kept
    write a file to to run treemix. 
    since treemix only accepts biallelic sites, also print out some other statistic numbers. For each position in the fasta sequences, if count of a base is less than Ncut, then it will be considered as a gap and not considered
    if output == None, write to 'output.treemix'
    if gapCut==1, then output only sites with no gaps.
    allStrand is True, means all strand in map file will be used
    threads is the number of CPUs to use when doing multiple processing
    '''
    # get mapfiles locations and their species
    print('read in input map files')
    if isinstance(maps, list):
        print('input maps is a list, use directly')
    elif isinstance(maps, str):
        print('input maps is a filename')
        maps = open(maps).readlines()
        maps = [e.split() for e in maps]
    else:
        print('input maps is not in correct format')
        return None
    print('there are ', len(maps), 'mapfiles')
    
    # get the relationship between mapfile and species
    ls_maps = [e[0] for e in maps]
    ls_mapOpen = [open(e) for e in ls_maps]
    ls_species = [e[1] for e in maps]
    print(len(set(ls_species)),'species in total')
    
    # prepare the output file
    if output is None:
        output = 'output.treemix'
    fout = open(output,'w')
    
    # get max allowed gap per site
    max_gap = _getMaxAllowedGaps(ls_maps,  allStrand=allStrand, gapCut=gapCut)
    print('the number of gaps per site should be smaller than',max_gap)
    
    # process siteKeep
    global sitesUse
    if siteKeep is None:
        sitesUse = False
        print('siteKeep not provided, use all sites')
    else:
        sitesUse = open(siteKeep).read()
        sitesUse = sitesUse.split()
        sitesUse = list(map(int,sitesUse))
        sitesUse = set(sitesUse)
        print('total sites in siteKeep is',len(sitesUse))
    
    
    
    for lines in zip(*ls_mapOpen):
        pass
    
    seqs = list(SeqIO.parse(file_fasta,'fasta'))
    print('total numbers of sequences in file_fasta', len(seqs))
    seqs = [e for e in seqs if e.id in dc_sample2group]
    print('number of sequences exist in file_groups that will be used',len(seqs))
    seq_ids = [e.id for e in seqs]
    seqs = [str(e.seq) for e in seqs]
    dc_seqs = dict(zip(seq_ids,seqs))
    
    seqlen = len(seqs[0])
    seqNum = len(seqs)
    print('there are', seqNum, 'sequences with the length of', seqlen, 'will be used. \n update the original group2sample, sample2group dictionary')
    dc_sample2group = {k:v for k,v in dc_sample2group.items() if k in seq_ids}
    dc_group2sample = {}
    for k,v in dc_sample2group.items():
        if v not in dc_group2sample:
            dc_group2sample[v] = []
        dc_group2sample[v].append(k)
    min_group_samples = min(len(v) for v in dc_group2sample.values())
    print('the smallest group have', min_group_samples, 'samples')
    if Ncut >= min_group_samples:
        print('\n\nwarning, Ncut is greater than the minimum sample count for a group, which means that the group might have no sequences left after filtering. Consider to change Ncut\n\n')
    
    ls_baseCounts = [Counter(e) for e in zip(*seqs)]
    print('first filter Ncut, change counts of ATCG to 0 if it is less than', Ncut)
    sites_changed = 0
    for e in ls_baseCounts:
        changed = False
        for k in 'ATCG':
            if e[k] < Ncut:
                e[k] = 0
                changed = True
        if changed:
            sites_changed += 1
    print(sites_changed, 'sites changed after filter Ncut')
    
    min_nongap = np.ceil(seqNum * gapCut)
    max_gap = seqNum - min_nongap
    print('filter gapCut, totally', seqNum, 'sequences, the allowed max gap per site is',max_gap)
    sites_keep = []
    for n, baseCounts in enumerate(ls_baseCounts):
        nongapCounts = sum(v for k,v in baseCounts.items() if k in 'ATCG')
        if nongapCounts >= min_nongap:
            sites_keep.append(n)
    print((seqlen - len(sites_keep))/seqlen, 'of sites removed after filter gapCut')
    
    print('filter  biallelic sites')
    sites_keep2 = []
    for n in sites_keep:
        baseCounts = ls_baseCounts[n]
        biallelic_bases = [e for e in 'ATCG' if baseCounts[e]>0]
        if len(biallelic_bases) == 2:
            sites_keep2.append((n,biallelic_bases))
    print(len(sites_keep2),'sites left as input for treemix')
    
    dc_group2sampleSeqs = {}
    for k in dc_group2sample:
        if k not in dc_group2sampleSeqs:
            dc_group2sampleSeqs[k] = []
        for v in dc_group2sample[k]:
            dc_group2sampleSeqs[k].append(dc_seqs[v])
    
    dc_group2sampleBiallic ={}
    dc_sites_keep2 = {e[0]:e[1] for e in sites_keep2}
    for k in dc_group2sampleSeqs:
        if k not in dc_group2sampleBiallic:
            dc_group2sampleBiallic[k] = []
        for n, bases in enumerate(zip(*dc_group2sampleSeqs[k])):
            if n in dc_sites_keep2:
                baseCounts = Counter(bases)
                bases = dc_sites_keep2[n]
                bases_c = [baseCounts[e] for e in bases]
                if len(bases_c) != 2:
                    print('something wrong here')
                dc_group2sampleBiallic[k].append(str(bases_c[0])+','+str(bases_c[1]))
    
    if output is None:
        output = file_fasta + '.treemix'
    outpath = os.path.dirname(output)
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    fout = open(output,'w')
    ls_group2sampleBiallic = []
    for k,v in dc_group2sampleBiallic.items():
        ls_group2sampleBiallic.append([k]+v)
    for line in zip(*ls_group2sampleBiallic):
        fout.write('\t'.join(line)+'\n')
    fout.close()
    
    #gzip the file
    os.system('gzip -f '+output)
    print('done')


description = '''
    file_fasta is a fasta file. each fasta id is a sample name that belongs to a population group
    file_groups is a file of relationship between sample name and population name
    if a sample is not in file_groups, then do not output its result
    write a file to to run treemix. 
    since treemix only accepts biallelic sites, also print out some other statistic numbers. For each position in the fasta sequences, if count of a base is less than Ncut, then it will be considered as a gap and not considered
    if output == None, write to file_fasta+'.treemix'
    if gapCut==1, then output only sites with no gaps.
    '''
    

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input file of fasta sequences', required=True)
    parser.add_argument('-g','--groups',help = 'sample to group table', required = True)
    parser.add_argument('-o','--output',help = 'where to store the output file. the output file is a input for treemix',default=None)
    parser.add_argument('-N','--Ncut', help = 'minimum bases couts to use. if NCut = 4, (A3T5G5) will be changed to (A0T5G5), and since A is 0, the site is still consider biallelic. default 0', default = 0, type=int)
    parser.add_argument('-G','--gapCut', help = 'minimum non-gap base proportion at each site. default 1, which means no gap is allowed.', default = 1, type=float)
    f = parser.parse_args()
    fasta2treemixInput(file_fasta=f.input, file_groups=f.groups, Ncut = f.Ncut, gapCut=f.gapCut, output=f.output)

