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
        #need to change accordingly here
        gene = items[-1].split('-')[0].replace('ID=', '')
        if scaf not in coding:
            coding[scaf] = {}


        i, j = list(map(int, items[3:5]))
        for k in range(i, j+1):
            coding[scaf][k] = gene


#2. filter header for coding regions

codingI = {}
with open(fhead) as fp:
    for i, line in enumerate(fp):
        items = line.strip().split()
        scaf, index = items[:2]
        index = int(index)
        try:
            cDict = coding[scaf]
        except KeyError:
            continue

        try:
            gene = cDict[index]
        except KeyError:
            continue

        try:
            codingI[gene].append(i)
        except KeyError:
            codingI[gene] = [i]


outdir = 'coding_fasta'
cmn.run('rm -r ' + outdir)
cmn.mkdir(outdir)


with open(ffa) as fp:
    for line in fp:
        line = line.strip()
        if line[0] == '>':
            defline = line[1:]
        else:
            seq = line
            for gene in codingI:
                dn = '%s/%s.fa' % (outdir, gene)
                geneSeq = ''.join([seq[i] for i in codingI[gene]])
                fasta = '>%s\n%s\n' % (defline, geneSeq)
                cmn.append_file(fasta, dn)

