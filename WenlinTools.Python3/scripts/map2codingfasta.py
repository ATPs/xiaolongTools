import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


fmap = sys.argv[1]
fhead = sys.argv[2]
fgff = sys.argv[3]

sp = cmn.lastName(fmap).split('_')[0]

#1. read in gff definition

coding = {}

with open(fgff) as fp:
    for line in fp:
        if line.strip() == '':
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

codingI = set([])
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
            codingI.add(i)



#3. filtering map for fasta
seqs = [[], []]

with open(fmap) as fp:
    for i, line in enumerate(fp):
        if i in codingI:
            a, b = line.strip().split()
            seqs[0].append(a)
            seqs[1].append(b)


dn = '%s_coding.fa' % sp
fasta = ['>%s_cp%s\n%s\n' % (sp, i+1, ''.join(seq)) for i, seq in enumerate(seqs)]
cmn.write_file(''.join(fasta), dn)

