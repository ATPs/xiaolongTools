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


rdict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        '-': 'N',
        '*': 'N'
        }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()


    lines = cmn.file2lines(fn)
    if lines[0][0] == '>':
        defline = lines[0]
        seq = ''.join(lines[1:])
    else:
        defline = '>' + cmn.lastName(fn)
        seq = ''.join(lines)

    seqlist = []
    for i in seq[::-1]:
        try:
            char = rdict[i]
        except KeyError:
            char = 'N'
        seqlist.append(char)
    seq = ''.join(seqlist)

    fasta = '%s_reverse\n%s\n' % (defline, seq)

    print(fasta)



