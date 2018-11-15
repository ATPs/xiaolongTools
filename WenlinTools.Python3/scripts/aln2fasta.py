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

def read_aln(fn):
    seqs = {}
    for line in cmn.file2lines(fn):
        name, seq = line.strip().split()
        seqs[name] = seq

    return seqs



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    seqDict = read_aln(fn)

    new = []
    for name in seqDict:
        seq = seqDict[name]
        new.append('>%s\n%s\n' % (name, seq))

    dn = cmn.lastName(fn) + '.fa'
    cmn.write_file(''.join(new), dn)


