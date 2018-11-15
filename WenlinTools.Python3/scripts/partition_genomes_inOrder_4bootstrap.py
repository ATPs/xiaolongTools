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
import random


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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        part = int(sys.argv[2])
    except:
        print("Usage: *.py *.fa Npart", file=sys.stderr)
        sys.exit()


    seqDict, length = read_fa(fn)

    positions = list(range(length))

    #random.shuffle(positions)

    bundle = length / part + 1

    label = cmn.lastName(fn).replace('.fa', '')
    outdir = 'BTsplitInOrder_%s' % label
    cmn.mkdir(outdir)

    for i in range(part):
        p = positions[i*bundle: ((i+1) * bundle)]

        dn = '%s/%s_%s.fa' % (outdir, label, i)
        with open(dn, 'w') as dp:
            for name in seqDict:
                seq = seqDict[name]
                newSeq = [seq[j] for j in p]
                dp.write('>%s\n' % name)
                dp.write(''.join(newSeq))
                dp.write('\n')



