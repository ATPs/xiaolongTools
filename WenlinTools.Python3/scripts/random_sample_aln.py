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
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, len(seq)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    N = 1000000 #get this much positions
    label = '%sM' % (N / 1000000)

    seqDict, length = read_fa(fn)

    if length < N:
        print('sequence length is shorter than 10K, exist!')
        sys.exit()

    positions = random.sample(list(range(length)), N)
    positions.sort()

    new = []
    for name in seqDict:
        seq = seqDict[name]

        newSeq = [seq[i] for i in positions]

        fasta = '>%s\n%s\n' % (name, ''.join(newSeq))
        new.append(fasta)

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '_rd%s.fa' % (label)
    cmn.write_file(''.join(new), dn)




