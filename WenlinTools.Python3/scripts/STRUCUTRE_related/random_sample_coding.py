#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import random

all_seq = {}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, len(seq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_rep():
    dn = 'rep.dict.pkl'
    if cmn.filexist(dn):
        print('loading repeats using precomputed data...')
        return cmn.pickle_read(dn)

    freps = cmn.cmd2lines('ls annotation_repeats/*.gff3')

    repdict = {}
    for frep in freps:
        for line in cmn.file2lines(frep):
    #if True:
    #    for line in cmn.cmd2lines('cat annotation_repeats/*.gff3| grep scaffold575_cov95'):
            items = line.strip().split()
            scaf = items[0]
            if scaf not in repdict:
                repdict[scaf] = set([])

            i, j = list(map(int, items[3:5]))
            repdict[scaf] = repdict[scaf] | set(range(i, j))
    cmn.pickle_write(repdict, dn)
    return repdict



def sample_regions(length, window, times):
    rlist = []
    while(len(rlist) <= times):
        i = random.randint(0, length-window)
        if len(rlist) == 0:
            rlist.append(i)
        else:
            if all([abs(i-each) > window for each in rlist]):
                rlist.append(i)
    genes = [(i, i+window) for i in rlist]
    return genes


if __name__=='__main__':

    ### random sample for structure run
    window = 10000 #sampled window
    #times = 10 #sample times
    times = int(sys.argv[1])

    seqDict, length = read_fa('all_coding_gaps_20_removed.ali')

    sampled_windows = sample_regions(length, window, times)

    ###1. read in gene regions
    #fn = '../introgression/pse.gff'
    #genes = []
    #count = 0
    #for line in cmn.file2lines(fn):
    #    items = line.strip().split()
    #    scaf = items[0]
    #    i, j = map(int, items[3:5])

    #    if j - i > window:
    #        count += 1
    #        genes.append((scaf, i,j))
    #print 'number of genes larger than %s: %s' % (window, count)

    ###3. prepare rep regions
    #repdict = read_rep()
    #repdict = {}


    ###4. take fasta
    final = {}
    for gene in sampled_windows:
        i, j = gene
        #fa = '../introgression/0_process_scaf/scaf2_fastas/%s.fa' % scaf
        #seqDict = read_fa(fa)
        #print 'parsing %s (%s, %s)' % (fa, i, j)
        #try:
        #    exclusion = repdict[scaf]
        #except:
        #    exclusion = set([])

        for name in seqDict:
            #seq = [char for index, char in enumerate(seqDict[name])
            #        if index not in exclusion and (i <= index <= j)]
            seq = seqDict[name][i:j]
            #seq = [char for index, char in enumerate(seq)
            #        if (index+i) not in exclusion]
            try:
                final[name] += seq
            except:
                final[name] = seq

    dn = 'sampled_seq_coding_t%s.fa' % times
    fasta = ['>%s\n%s\n' % (name, ''.join(final[name]))
            for name in final]
    cmn.write_file(''.join(fasta), dn)





