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
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

reverse_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        }

def reverse_strand(read):
    new = []
    for char in read[::-1]:
        try:
            a = reverse_dict[char]
        except KeyError:
            a = 'N'
        new.append(a)

    return ''.join(new)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fns= [os.path.abspath(each) for each in sys.argv[3:]]
        #KmerCut = int(sys.argv[1])
        KmerSize = int(sys.argv[1])
        Ncpu = int(sys.argv[2])
    except:
        print("Usage: *.py KmerCut KmerSize(19) Ncpu R1.fq R2.fq", file=sys.stderr)
        sys.exit()



    outlabel = cmn.lastName(fns[0]).split('_')[0]
    tmpDir = '%s_jf' % outlabel
    cmn.mkdir(tmpDir)
    #step1, run Jellyfish
    print('running Jellyfish to get Kmer count...')
    os.chdir(tmpDir)
    cmd = 'jellyfish count -m %s -t %s -s 10000000000 -c 8 --timing=jf.err --canonical ' % (KmerSize, Ncpu)
    cmd += ' '.join(fns)
    cmn.run(cmd)

    cmd = 'jellyfish histo mer_counts.jf > %smer_histo.txt' % KmerSize
    cmn.run(cmd)

    cmd = 'jellyfish dump -c mer_counts.jf > %smer_counts' % KmerSize
    cmn.run(cmd)

    #step2, filter out reads with high Kmers
    #for each read to be filtered out, require the high coverage Kmer appear more than once
    #it very high Kmer in a read, dump it directly

    #read in
    veryHighSet = set([])
    HighSet = set([])
    cut1 = 10000
    cut2 = 1000
    with open('%smer_counts' % KmerSize) as fp:
        for line in fp:
            seq, count = line.strip().split()
            count = int(count)
            if count > cut1:
                veryHighSet.add(seq)
            else:
                if count > cut2:
                    HighSet.add(seq)

    os.chdir('..')
    #parsing
    for fn in fns:
        label = cmn.lastName(fn).replace('.fastq', '').replace('.fq', '')
        dstat = '%s_Kmer.stat' % label
        with open(dstat, 'w') as dpStat, open(fn) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 0:
                    ID = line[1:].strip()
                    count1 = 0
                    count2 = 0
                elif i % 4 == 1:
                    seq = line.strip()
                    rseq = reverse_strand(seq)
                    for i in range(len(seq) - KmerSize):
                        segment1 = seq[i:i+KmerSize]
                        segment2 = rseq[i:i+KmerSize]
                        for segment in [segment1, segment2]:
                            if segment in veryHighSet:
                                count1 += 1
                            else:
                                count2 += 1
                    if count1 != 0 or count2 != 0:
                        line = '%s\t%s\t%s\t%s\n' % (ID, count1, count2, seq)
                        dpStat.write(line)



    #TODO: use quake to fix those Kmers only appear once


