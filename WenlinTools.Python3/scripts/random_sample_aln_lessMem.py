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
        length = int(sys.argv[2])
    except:
        print("Usage: *.py fa alignment_length", file=sys.stderr)
        sys.exit()

    N = 10000000 #get this much positions
    label = '%sM' % (N / 1000000)

    #seqDict, length = read_fa(fn)

    #if length < N:
    #    print 'sequence length is shorter than 30K, exist!'
    #    sys.exit()

    positions = random.sample(list(range(length)), N)

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '_rd%s.fa' % (label)
    with open(dn, 'w') as dp:
        with open(fn) as fp:
            for line in fp:
                if line[0] == '>':
                    dp.write(line)
                else:
                    seq = line.strip()
                    newSeq = [seq[i] for i in positions]
                    seq = ''.join(newSeq) + '\n'
                    dp.write(seq)




