import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_exo(fn):
    #scaffold78778_cov57	exonerate:protein2genome:local	cds	35808	35945	.	-	.	HMEL015204-RA
    adict = {}
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        i, j = list(map(int, items[3:5]))
        scaf = items[0]
        protein = items[-1]
        if scaf not in adict:
            adict[scaf] = {}

        try:
            adict[scaf][protein].append((i,j))
        except KeyError:
            adict[scaf][protein] = [(i, j)]
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
maplist = cmn.file2lines(sys.argv[1])
#fmap = sys.argv[1]
fhead = sys.argv[2]
#fgff = sys.argv[3]
fexo = sys.argv[3]

#sp = cmn.lastName(fmap).split('_')[0]

#1. read in gff definition

#coding = cmn.pickle_read(fgff)
coding_dict = read_exo(fexo)
scafs = set(coding_dict.keys())

#2. filter header for coding regions

codingI = {}
with open(fhead) as fp:
    for i, line in enumerate(fp):
        items = line.strip().split()
        scaf, index = items[:2]
        index = int(index)
        try:
            cIndex_dict = coding_dict[scaf]
        except KeyError:
            continue

        for protein in cIndex_dict:
            alist = cIndex_dict[protein]
            if any([i<=index<=j for i,j in alist]):
                codingI[i] = protein



#3. filtering map for fasta
#seqs = [[], []]
outdir = 'gene_fasta'
cmn.mkdir(outdir)


for fmap in maplist:
    sp = cmn.lastName(fmap).split('_')[0]

    prot_dict = {}
    with open(fmap) as fp:
        for i, line in enumerate(fp):
            try:
                prot = codingI[i]
            except KeyError:
                continue

            #reach here if protein found
            a, b = line.strip().split()
            try:
                prot_dict[prot].append((a, b))
            except KeyError:
                prot_dict[prot] = [(a, b)]

    for prot in prot_dict:
        dn = '%s/%s.fa' % (outdir, prot)
        alist = prot_dict[prot]
        new = []
        for i in range(2):
            seq = ''.join([each[i] for each in alist])
            name = '%s_cp%s' % (sp, i+1)
            fasta = '>%s\n%s\n' % (name, seq)
            cmn.append_file(fasta, dn)

    print('finish processing %s' % fmap)
