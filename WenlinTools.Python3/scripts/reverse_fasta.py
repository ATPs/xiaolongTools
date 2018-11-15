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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split()
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq.upper()
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py aln.fa", file=sys.stderr)
        sys.exit()


    seqDict = read_fa(fn)

    newDict = {key: ''.join([rdict[char] for char in seqDict[key][::-1]])
            for key in seqDict}

    dn = cmn.lastName(fn) + '.reverse'
    fastas = ['>%s_reverse\n%s\n' % (name, newDict[name])
            for name in newDict]
    cmn.write_file(''.join(fastas), dn)


