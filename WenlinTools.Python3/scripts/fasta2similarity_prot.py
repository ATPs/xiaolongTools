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
gapChars = set(['X', '-', '*'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    alist = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        alist.append(defline)
        adict[defline] = seq
    return adict, alist



def compute_identity(a, b):
    #print a
    #print b
    seqlength = len(a)
    #if len(b) < len(a):
    #    seqlength = len(b)
    if len(a) != len(b):
        return 0, 0, 0

    check_list = [a[i] == b[i] for i in range(seqlength)
            if a[i] not in gapChars and b[i] not in gapChars]
    if len(check_list) == 0:
        return 0, 0, 0
    Nsame = sum(check_list)
    identity = float(Nsame) / len(check_list)
    return identity, Nsame, len(check_list)


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()

    #snice this code was used to get barcode for Nick, add a function of checking alignment length

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '.clw'

    seqDict, orderlist = read_fa(fn)

    keys = list(seqDict.keys())
    for i, name1 in enumerate(keys):
        seq1 = seqDict[name1]
        for name2 in keys[i+1:]:
            seq2 = seqDict[name2]
            idt, Nsame, countN = compute_identity(seq1, seq2)
            print(name1, name2, countN, (countN - Nsame), idt)

