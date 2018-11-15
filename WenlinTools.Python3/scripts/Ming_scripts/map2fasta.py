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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.map", file=sys.stderr)
        sys.exit()


    label = cmn.lastName(fn).split('_')[0]
    dn = cmn.lastName(fn).replace('.map', '') + '_m2s.fa'

    print('parsing for %s' % label)
    print('would save result to %s' % dn)

    a = []
    b = []

    with open(fn) as fp:
        for line in fp:
            i, j = line.strip().split()
            a.append(i)
            b.append(j)


    fasta = '>%s_cp1\n%s\n>%s_cp2\n%s\n' % (label, ''.join(a), label, ''.join(b))
    cmn.write_file(fasta, dn)



