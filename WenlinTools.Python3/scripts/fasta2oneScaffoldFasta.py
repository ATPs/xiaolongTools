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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fhead, taken_scaf = sys.argv[1:]
    except:
        print("Usage: *.py fa fhead scafname", file=sys.stderr)
        sys.exit()


    #fhead = '/work/biophysics/mtang/SNP_calling/indexed_references/Junonia_v2_scaf.header'
    #fhead = 'Calycopis_cecrops_assembly_V1.1_scaf.header'
    print('loading header info...')
    headDict = {}
    with open(fhead) as fp:
        for i, line in enumerate(fp):
            scaf, index = line.strip().split()
            if scaf != taken_scaf:
                continue
            try:
                headDict[scaf].append(i)
            except KeyError:
                headDict[scaf] = [i]

    print('finish loading header, begin parsing fasta...')
    #outdir = '%s_scafs' % cmn.lastName(fn)
    #cmn.mkdir(outdir)

    seqDict = read_fa(fn)
    for scaf in headDict:
        indexes = headDict[scaf]
        new = []
        for name in seqDict:
            seq = seqDict[name]
            newSeq = ''.join([seq[i] for i in indexes])
            fasta = '>%s\n%s\n' % (name, newSeq)
            new.append(fasta)

        dn = '%s_%s.fa' % (cmn.lastName(fn), taken_scaf)
        cmn.write_file(''.join(new), dn)




