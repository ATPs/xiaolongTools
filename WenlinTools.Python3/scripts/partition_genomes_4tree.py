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
#import random


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    with open(fa) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '>':
                defline = line[1:]
            else:
                seq = line
                try:
                    adict[defline].append(seq)
                except KeyError:
                    adict[defline] = [seq]
    keys = list(adict.keys())

    for key in keys:
        seq = ''.join(adict[key])
        adict[key] = seq


    return adict, len(seq)


def parse_header_file(fn):
    mapDict = {}
    lengthDict = {}
    with open(fn) as fp:
        for i, line in enumerate(fp):
            scaf, index = line.strip().split()
            if 'mito' in scaf:
                continue
            try:
                mapDict[scaf].append(i)
            except KeyError:
                mapDict[scaf] = [i]
            try:
                lengthDict[scaf] += 1
            except KeyError:
                lengthDict[scaf] = 1
    return lengthDict, mapDict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gapChars = set(list('-N'))


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        fhead = sys.argv[2]
        lengthCut = int(sys.argv[3])
    except:
        print("Usage: *.py ass.fa *.header 5000", file=sys.stderr)
        print("take the scaffolds larger than [5000] and partion them into pieces of [5000]", file=sys.stderr)
        sys.exit()

    scafLength, mapdict = parse_header_file(fhead)
    #sorted_scafs = sorted(scafLength.keys(), key=lambda x:scafLength[x])
    filtered_scafs = [scaf for scaf in scafLength
            if scafLength[scaf] >= lengthCut]

    seqDict, length = read_fa(fn)
    depth = float(len(seqDict))

    outdir = 'split4treeWindow%s' % lengthCut

    cmn.mkdir(outdir)
    for scaf in filtered_scafs:
        indexes = mapdict[scaf]
        scaf_length = scafLength[scaf]
        i = 0
        while(i < scaf_length):
            j = i + lengthCut
            dn = '%s/%s_%s-%s.fa' % (outdir, scaf, i, j)
            subPs = indexes[i:j]
            with open(dn, 'w') as dp:
                for name in seqDict:
                    seq = seqDict[name]
                    newSeq = [seq[k] for k in subPs]
                    dp.write('>%s\n' % name)
                    dp.write(''.join(newSeq))
                    dp.write('\n')
            i += lengthCut


