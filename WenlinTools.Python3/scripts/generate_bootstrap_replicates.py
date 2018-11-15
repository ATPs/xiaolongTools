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
        N = int(sys.argv[2])
    except:
        print("Usage: *.py faln 100", file=sys.stderr)
        sys.exit()

    dnlabel = cmn.lastName(fn)
    seqDict, length = read_fa(fn)

    for rep in range(N):
        dn = '%s_rep%s' % (dnlabel, rep)
        print('processing %s' % dn)
        taken_indexes = []
        for _ in range(length):
            taken_indexes.append(random.randint(0, length - 1))

        with open(dn, 'w') as dp:
            for name in seqDict:
                seq = seqDict[name]
                newSeq = [seq[i] for i in taken_indexes]
                fa = '>%s\n%s\n' % (name, ''.join(newSeq))
                dp.write(fa)


