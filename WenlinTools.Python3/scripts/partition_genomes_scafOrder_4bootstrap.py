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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        fhead = sys.argv[2]
        part = int(sys.argv[3])
    except:
        print("Usage: *.py ass.fa *.header 100", file=sys.stderr)
        sys.exit()

    scafLength, mapdict = parse_header_file(fhead)
    sorted_scafs = sorted(list(scafLength.keys()), key=lambda x:scafLength[x])

    seqDict, length = read_fa(fn)

    #random.shuffle(positions)

    bundle = length / part + 1

    label = cmn.lastName(fn).replace('.fa', '')
    outdir = 'BTsplitScafOrder_%s' % label
    cmn.mkdir(outdir)

    positions = []
    count = 0
    for scaf in sorted_scafs:
        #check for each bundle
        indexes = mapdict[scaf]
        for index in indexes:
            positions.append(index)
            if len(positions) >= bundle:
                dn = '%s/%s_%s.fa' % (outdir, label, count)
                with open(dn, 'w') as dp:
                    for name in seqDict:
                        seq = seqDict[name]
                        newSeq = [seq[j] for j in positions]
                        dp.write('>%s\n' % name)
                        dp.write(''.join(newSeq))
                        dp.write('\n')
                count += 1
                positions = []

    if positions != []:
        if True:
            if True:
                dn = '%s/%s_%s.fa' % (outdir, label, count)
                with open(dn, 'w') as dp:
                    for name in seqDict:
                        seq = seqDict[name]
                        newSeq = [seq[j] for j in positions]
                        dp.write('>%s\n' % name)
                        dp.write(''.join(newSeq))
                        dp.write('\n')





