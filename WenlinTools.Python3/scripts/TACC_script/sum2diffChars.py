#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.clw", file=sys.stderr)
        sys.exit()

    adict = {}
    with open(fn) as fp:
        for line in fp:
            exon, sp, seq = line.strip().split()
            sp = sp.split('.')[0]
            exon_key = exon[:-2]
            if sp not in adict:
                adict[sp] = {}
            try:
                adict[sp][exon_key].append(seq)
            except KeyError:
                adict[sp][exon_key] = [seq]
    new = []
    sumDict = {}
    for sp in adict:
        for exon in adict[sp]:
            seqs = adict[sp][exon]
            N = sum([seqs[0][i] != seqs[1][i]
                for i in range(len(seqs[0]))])
            new.append('%s\t%s\t%s\n' % (sp, exon, N))
            try:
                sumDict[sp] += N
            except KeyError:
                sumDict[sp] = N

    for sp in sumDict:
        new.append('spSum\t%s\t%s\n' % (sp, sumDict[sp]))
    dn = cmn.lastName(fn).replace('.sum', '') + '.diff'
    cmn.write_file(''.join(new), dn)


