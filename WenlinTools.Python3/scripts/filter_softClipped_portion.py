#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pysam
import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def cut_ends(reflength, aligned):
    #i is read, j is ref
    minN, maxN = [], []
    readLength = len(aligned)

    #print aligned
    for i, j in aligned:
        if i != None and j != None:
            maxN = [i, j]

            if minN == []:
                minN = [i, j]

    leftN = 0
    rightN = readLength
    if minN[1] - minN[0] < 0:
        leftN = minN[0] - minN[1]

    if maxN[1] + (readLength - maxN[0]) > reflength:
        rightN = maxN[0] + (reflength - maxN[1])

    #if rightN != readLength:
    #    print 'cutting right: %s' % str(aligned)
    #    print aligned[leftN: rightN]
    return aligned[leftN: rightN]



def check_aligned_range(aligned):
    i = 0
    while(i < len(aligned)):
        a, b = aligned[i]
        if a != None and b != None:
            break
        i += 1

    j = len(aligned) - 1
    while(j > 0):
        a, b = aligned[j]
        if a != None and b != None:
            break
        j -= 1

    j += 1
    return i, j


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.sam", file=sys.stderr)
        sys.exit()


    samfile = pysam.AlignmentFile(fn)
    reflength = samfile.lengths[0]

    dnlabel = cmn.lastName(fn)
    dn = 'noSC_%s' % dnlabel
    cmd = 'grep ^@ %s > %s' % (fn, dn)
    cmn.run(cmd)

    Nbad = 0
    Ntotal = 0
    goods = []
    cut_ranges = {}
    for i, record in enumerate(samfile):
        if record.is_secondary or record.is_unmapped:
            continue

        Ntotal += 1
        aligned = record.aligned_pairs
        if len(aligned) == 0:
            continue

        i, j = check_aligned_range(aligned)
        q = record.qual
        record.seq = record.seq[i:j]
        record.qual = q[i:j]

        goods.append(record)

    #print 'filter out %s bad reads out of total %s' % (Nbad, Ntotal)

    dn2 = 'tmp_work_%s' % dnlabel
    with pysam.AlignmentFile(dn2, 'w', header=samfile.header) as outf:
        for record in goods:
            outf.write(record)

    cmn.run('cat %s >> %s; rm %s' % (dn2, dn, dn2))


