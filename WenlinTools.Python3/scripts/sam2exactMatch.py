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
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
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


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fref = sys.argv[1:]
    except:
        print("Usage: *.py sam ref.fa", file=sys.stderr)
        sys.exit()

    seqDict = read_fa(fref)

    dnH = cmn.lastName(fn).replace('.sam', '') + '.matchStat'

    hdict = {}
    samfile = pysam.AlignmentFile(fn)
    for record in samfile:
        if record.is_unmapped:
            continue
        if record.is_secondary:
            continue

        scaf = record.reference_name
        aligns = record.get_aligned_pairs()
        query = record.query_sequence
        #sbjct = record.reference_sequence
        sbjct = seqDict[scaf]
        Nmismatch = sum([query[queryI] != sbjct[sbjctI]
            for queryI, sbjctI in aligns
            if queryI != None and sbjctI != None])
        N = len([each for each in aligns
            if None not in each])

        if Nmismatch > 1:
            continue

        if N >= 40 and float(N) / len(aligns) > 0.8:

            if scaf not in hdict:
                hdict[scaf] = [0, 0]
            hdict[scaf][0] += 1
            hdict[scaf][1] += N


    alist = sorted(hdict, key=lambda x: hdict[x])
    with open(dnH, 'w') as dp:
        for key in alist:
            line = '%s\t%s\t%s\n' % (key, hdict[key][0], hdict[key][1])
            dp.write(line)






