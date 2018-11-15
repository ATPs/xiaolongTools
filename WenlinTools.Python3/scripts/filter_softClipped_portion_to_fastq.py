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

reverse_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        }

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


def record2name(record, showRef=False):
    Oname = record.query_name
    if record.is_read1:
        name = Oname + '/1'
        i = 0
    elif record.is_read2:
        name = Oname + '/2'
        i = 1
    if showRef:
        name = '%s|%s' % (name, record.reference_name)
    return Oname, i, name



def shrink_ends(i, j, length):
    N = 6
    i = i + N
    if i > length:
        print('error here! %s larger than length %s' % (i, length))
        i = length

    j = j - N
    if j < 0:
        print('error here! %s smaller than 0' % j)
        j = 0
    return i, j


def reverse_strand(read):
    new = []
    for char in read[::-1]:
        try:
            a = reverse_dict[char]
        except KeyError:
            a = 'N'
        new.append(a)

    return ''.join(new)


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        outlabel = sys.argv[2]
    except:
        print("Usage: *.py *.sam outlabel", file=sys.stderr)
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
    rdict = {}
    for i, record in enumerate(samfile):
        if record.is_secondary or record.is_unmapped:
            continue

        Ntotal += 1
        aligned = record.aligned_pairs
        if len(aligned) == 0:
            continue

        i, j = check_aligned_range(aligned)
        if j - i < len(record.seq):
            #i, j = shrink_ends(i, j, len(aligned))
            Nbad += 1

        name, iii, subName = record2name(record)
        q = record.qual[i:j]
        s = record.seq[i:j]
        if iii == 1:
            s = reverse_strand(s)
            q = q[::-1]
        fq = '@%s\n%s\n+\n%s\n' % (subName, s, q)
        print(name, subName, i, j, record.seq, s)
        if name not in rdict:
            rdict[name] = [None, None]

        rdict[name][iii] = fq

    print('Nbad: %s; Ntotal %s;' % (Nbad, Ntotal))
    fq1 = '%s_R1.fq' % outlabel
    fq2 = '%s_R2.fq' % outlabel
    fsingle = '%s_singleton.fq' % outlabel
    for fn in [fq1,fq2,fsingle]:
        cmn.run('rm %s' % fn)

    for name in rdict:
        alist = rdict[name]
        if alist.count(None) == 0:
            cmn.append_file(alist[0], fq1)
            cmn.append_file(alist[1], fq2)
        else:
            for each in alist:
                if each != None:
                    cmn.append_file(each, fsingle)

