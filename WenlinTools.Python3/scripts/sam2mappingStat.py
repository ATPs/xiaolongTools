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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py sam", file=sys.stderr)
        sys.exit()

    dnH = cmn.lastName(fn).replace('.sam', '') + '.HighQmapStat'
    dnL = cmn.lastName(fn).replace('.sam', '') + '.mapStat'

    rdict = {}
    hdict = {}
    samfile = pysam.AlignmentFile(fn)
    for record in samfile:
        if record.is_unmapped:
            continue

        scaf = record.reference_name
        aligns = record.get_aligned_pairs()
        N = len([each for each in aligns
            if None not in each])

        if scaf not in rdict:
            rdict[scaf] = [0, 0]
        rdict[scaf][0] += 1
        rdict[scaf][1] += N

        if record.is_secondary:
            continue

        if N >= 40 and float(N) / len(aligns) > 0.8:

            if scaf not in hdict:
                hdict[scaf] = [0, 0]
            hdict[scaf][0] += 1
            hdict[scaf][1] += N

    alist = sorted(rdict, key=lambda x: rdict[x][1])
    with open(dnL, 'w') as dp:
        for key in alist:
            line = '%s\t%s\t%s\n' % (key, rdict[key][0], rdict[key][1])
            dp.write(line)


    alist = sorted(hdict, key=lambda x: hdict[x])
    with open(dnH, 'w') as dp:
        for key in alist:
            line = '%s\t%s\t%s\n' % (key, hdict[key][0], hdict[key][1])
            dp.write(line)






