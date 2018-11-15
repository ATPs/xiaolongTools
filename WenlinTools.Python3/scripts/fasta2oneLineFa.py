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
    alist = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        alist.append(defline)
        adict[defline] = seq
    return adict, alist




if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()

    #snice this code was used to get barcode for Nick, add a function of checking alignment length

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '_1line.fa'

    seqDict, orderlist = read_fa(fn)

    with open(dn, 'w') as dp:
        for name in orderlist:
            seq = seqDict[name]
            dp.write('>%s\n%s\n' % (name, seq))
