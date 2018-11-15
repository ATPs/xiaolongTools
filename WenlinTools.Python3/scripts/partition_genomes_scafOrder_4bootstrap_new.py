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
        part = int(sys.argv[3])
        gapCut = float(sys.argv[4])
    except:
        print("Usage: *.py ass.fa *.header 100 0.3", file=sys.stderr)
        sys.exit()

    scafLength, mapdict = parse_header_file(fhead)
    sorted_scafs = sorted(list(scafLength.keys()), key=lambda x:scafLength[x])

    seqDict, length = read_fa(fn)
    depth = float(len(seqDict))

    #random.shuffle(positions)
    positions = []
    for scaf in sorted_scafs:
        indexes = mapdict[scaf]
        for index in indexes:
            gapN = len([1 for name in seqDict
                if seqDict[name][index] in gapChars])

            if gapN / depth < gapCut:
                positions.append(index)

    bundle = len(positions) / part + 1

    label = cmn.lastName(fn).replace('.fa', '')
    outdir = 'BTsplitScafOrderEq%sGap_%s' % (gapCut, label)
    cmn.mkdir(outdir)

    count = 0
    for i in range(part):
        subPs = positions[i*bundle : (i+1) * bundle]
        dn = '%s/%s_%s.fa' % (outdir, label, count)
        with open(dn, 'w') as dp:
            for name in seqDict:
                seq = seqDict[name]
                newSeq = [seq[j] for j in subPs]
                dp.write('>%s\n' % name)
                dp.write(''.join(newSeq))
                dp.write('\n')
        count += 1






