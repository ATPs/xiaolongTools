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
from collections import Counter


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
gaps = set(['X','N','-'])

def tell_informatic_positions(seqDict, length):
    p = []
    for i in range(length):
        chars = [seqDict[key][i] for key in seqDict]
        chars = [char for char in chars
                if char not in gaps]
        #TODO: currently, do not care about gap fraction here
        count_dict = Counter(chars)

        if all([count_dict[char] >= 4 for char in count_dict]):
            p.append(i)
    return p



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    seqDict, length = read_fa(fn)

    positions = tell_informatic_positions(seqDict, length)

    new = []
    for name in seqDict:
        seq = seqDict[name]

        newSeq = [seq[i] for i in positions]

        fasta = '>%s\n%s\n' % (name, ''.join(newSeq))
        new.append(fasta)

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '_infoP.fa'
    cmn.write_file(''.join(new), dn)




