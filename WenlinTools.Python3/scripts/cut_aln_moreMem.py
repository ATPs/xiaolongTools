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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    print(fastas)
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join([line.strip() for line in lines[1:]])
        #seq = seq.replace('N', '-')
        adict[defline] = seq
    return adict




if __name__=='__main__':
    #options=parse_options()
    try:
        fn, Range = sys.argv[1:3]
        i, j = list(map(int, Range.split('-')))
    except:
        print("Usage: *.py aln 0-10000", file=sys.stderr)
        sys.exit()

    new = []
    seqDict = read_fa(fn)
    for name in seqDict:
        seq = seqDict[name]
        fasta = '>%s\n%s\n' % (name, seq[i:j])
        new.append(fasta)

    dn = '%s_%s.fa' % (cmn.lastName(fn).replace('.fa', ''), Range)
    cmn.write_file(''.join(new), dn)



