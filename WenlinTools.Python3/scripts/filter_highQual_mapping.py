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
    dn = 'highQ_%s' % dnlabel
    cmd = 'grep ^@ %s > %s' % (fn, dn)
    cmn.run(cmd)

    Nbad = 0
    Ntotal = 0
    goods = []
    for record in samfile:
        if record.is_secondary:
            continue

        Ntotal += 1
        aligned = record.aligned_pairs
        if len(aligned) == 0:
            continue

        aligned = cut_ends(reflength, aligned)
        Naligned = sum([each[1] != None and each[0] != None for each in aligned])
        seq = record.query_sequence
        length = len(seq)
        name = record.query_name
        p = float(Naligned)/length
        #print name, p
        if p >= 0.8 and Naligned >= 40:
            goods.append(record)
        else:
            Nbad += 1

    print('filter out %s bad reads out of total %s' % (Nbad, Ntotal))

    dn2 = 'tmp_work_%s' % dnlabel
    with pysam.AlignmentFile(dn2, 'w', header=samfile.header) as outf:
        for record in goods:
            outf.write(record)

    cmn.run('cat %s >> %s; rm %s' % (dn2, dn, dn2))


