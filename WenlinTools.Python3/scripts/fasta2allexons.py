#!/usr/bin/env python

import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


ffa = sys.argv[1]
#fhead = sys.argv[2]
fgff = sys.argv[2]

sp = cmn.lastName(ffa).split('_')[0]

#1. read in gff definition

coding = {}
cds_count = {}

gene_count = {}
with open(fgff) as fp:
    for line in fp:
        if line.strip() == '':
            continue
        if line[0] == '#':
            continue

        items = line.strip().split()
        scaf = items[0].split()[0]

        if items[2] == 'gene':
            scafN = scaf.split('_')[0].replace('scaffold', '')
            if scafN not in gene_count:
                gene_count[scafN] = 1
            else:
                gene_count[scafN] += 1

        if items[2] != 'CDS':
            continue


        #print scaf
        gene = items[-1]
        if '=' in gene:
            gene = '%s.%s' % (scafN, gene_count[scafN])

        if gene not in cds_count:
            cds_count[gene] = 1
        else:
            cds_count[gene] += 1

        if scaf not in coding:
            coding[scaf] = {}


        i, j = list(map(int, items[3:5]))
        exon = '%s.e%s' % (gene, cds_count[gene])
        for k in range(i, j+1):
            coding[scaf][k] = exon


dn = '%s_exons.fa' % (cmn.lastName(ffa).replace('.fasta', '').replace('.fa', ''))

seqDict = {}
seq = []
defline = ''
with open(ffa) as fp:
    for line in fp:
        if line[0] == '>':
            if len(seq) != 0 and (defline in coding):
                seqDict[defline] = ''.join(seq)

            defline = line[1:].strip().split()[0]
            seq = []
        else:
            seq.append(line.strip())

seqDict[defline] = ''.join(seq)


with open(dn, 'w') as dp:
    for scaf in coding:
        iDict = coding[scaf]
        seq = seqDict[scaf]
        rdict = {}
        indexes = list(iDict.keys())
        indexes.sort()
        for index in indexes:
            exon = iDict[index]
            char = seq[index - 1]
            try:
                rdict[exon].append(char)
            except KeyError:
                rdict[exon] = [char]

        keys = list(rdict.keys())
        keys = sorted(keys, key=lambda x: int(x.split('e')[-1]))
        for exon in keys:
            fasta = '>%s\n%s\n' % (exon, ''.join(rdict[exon]))
            dp.write(fasta)

