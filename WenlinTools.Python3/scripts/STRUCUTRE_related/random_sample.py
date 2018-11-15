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
    global all_seq
    try:
        seqdict = all_seq[fa]
        return seqdict
    except KeyError:
        pass
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    all_seq[fa] = adict
    return adict

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



if __name__=='__main__':

    ### random sample for structure run
    window = 10000 #sampled window
    #times = 10 #sample times
    times = int(sys.argv[1])

    ###1. read in gene regions
    fn = '../introgression/pse.gff'
    genes = []
    count = 0
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        scaf = items[0]
        i, j = list(map(int, items[3:5]))

        if j - i > window:
            count += 1
            genes.append((scaf, i,j))
    print('number of genes larger than %s: %s' % (window, count))

    ###3. prepare rep regions
    repdict = read_rep()
    #repdict = {}

    ###2. random sample gene regions of 10K
    selected_genes = random.sample(genes, times)
    selected = []
    for gene in selected_genes:
        scaf, i, j = gene
        #try:
        #    exclusion = repdict[scaf]
        #except:
        #    exclusion = set([])

        randomI = random.randint(i, j-window)
        #while(randomI in exclusion):
        #    randomI = random.randint(i, j-window)

        randomJ = randomI + window
        selected.append((scaf, randomI, randomJ))

    ###4. take fasta
    final = {}
    for gene in selected:
        scaf, i, j = gene
        fa = '../introgression/0_process_scaf/scaf2_fastas/%s.fa' % scaf
        seqDict = read_fa(fa)
        print('parsing %s (%s, %s)' % (fa, i, j))
        try:
            exclusion = repdict[scaf]
        except:
            exclusion = set([])

        for name in seqDict:
            #seq = [char for index, char in enumerate(seqDict[name])
            #        if index not in exclusion and (i <= index <= j)]
            seq = seqDict[name][i:j]
            seq = [char for index, char in enumerate(seq)
                    if (index+i) not in exclusion]
            try:
                final[name] += seq
            except:
                final[name] = seq

    dn = 'sampled_seq_t%s.fa' % times
    fasta = ['>%s\n%s\n' % (name, ''.join(final[name]))
            for name in final]
    cmn.write_file(''.join(fasta), dn)





