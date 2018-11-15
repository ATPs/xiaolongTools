#!/usr/bin/env python

import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


ffa = sys.argv[1]
fhead = sys.argv[2]
fgff = sys.argv[3]

sp = cmn.lastName(ffa).split('_')[0]

#1. read in gff definition

coding = {}

with open(fgff) as fp:
    for line in fp:
        if line.strip() == '':
            continue
        if line[0] == '#':
            continue

        items = line.strip().split()
        if items[2] != 'CDS':
            continue

        scaf = items[0]
        if scaf not in coding:
            coding[scaf] = set([])

        i, j = list(map(int, items[3:5]))

        aset = set(range(i, j+1))

        coding[scaf] = coding[scaf] | aset


#2. filter header for coding regions

codingI = []
with open(fhead) as fp:
    for i, line in enumerate(fp):
        items = line.strip().split()
        scaf, index = items[:2]
        index = int(index)
        try:
            cIndex = coding[scaf]
        except KeyError:
            continue

        if index in cIndex:
            codingI.append(i)



#3. filtering map for fasta
seqs = [[], []]
dn = '%s_coding.fa' % sp
dp = open(dn, 'w')

with open(ffa) as fp:
    for line in fp:
        line = line.strip()
        if line[0] == '>':
            defline = line[1:]
            dp.write('>%s\n' % defline)
        else:
            coding_seq = [line[i] for i in codingI]
            coding_seq.append('\n')
            dp.write(''.join(coding_seq))
dp.close()

