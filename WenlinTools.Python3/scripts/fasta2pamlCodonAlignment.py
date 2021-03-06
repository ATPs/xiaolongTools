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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        #defline, seq = each.strip().split()
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:]).replace('N', '-')
        adict[defline] = seq
    return adict, len(seq)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def codon_seq(seq):
    seq = seq.replace('-', 'N')

    new = []
    for i in range(0, len(seq), 3):
        new.append(seq[i:i+3])

    return ' '.join(new)


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    seqDict, length = read_fa(fn)

    new = ['%s\t%s' % (len(seqDict), length)]
    for name in seqDict:
        new.append('%s        %s' % (name, codon_seq(seqDict[name])))


    dn = cmn.lastName(fn) + '.nuc'
    cmn.write_lines(new, dn)


