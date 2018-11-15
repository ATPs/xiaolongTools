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
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


def group_seqDict(adict):
    rdict = {}
    for name in adict:
        key = name.split('_')[0]
        seq = adict[name]
        try:
            rdict[key].append(seq)
        except:
            rdict[key] = [seq]
    return rdict


if __name__=='__main__':
    #fn = 'all_coding_gaps_20_removed.ali'
    try:
        fn = sys.argv[1]
    except:
        print('*.py fa', file=sys.stderr)
        sys.exit()

    seqDict = read_fa(fn)

    seqDict = group_seqDict(seqDict)

    GCset = set(list('GC'))
    #aset = set(cmn.getid('eubuleR1'))

    #spSeqs = [seqDict[name] for name in seqDict
    #        if name.split('_')[0] == '3314']

    #seqDict = {name: seqDict[name] for name in seqDict
    #        if name.split('_')[0] in aset}

    new = []
    for name in seqDict:
        seq1, seq2 = seqDict[name]
        length = len(seq1)
        line = [name]
        Ngap = seq1.count('-')
        Ngc = seq1.count('G') + seq1.count('C') + seq2.count('G') + seq2.count('C')
        Ngc = Ngc / 2

        print(name, Ngc, length, Ngap)

    #dn = 'hetero_check.out'
    #cmn.write_lines(new, dn)

