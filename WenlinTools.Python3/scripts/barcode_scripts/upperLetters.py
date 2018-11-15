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
    alist = []
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        alist.append(defline)
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, alist


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    seqDict, orderlist = read_fa(fn)

    new = ['>%s\n%s\n' % (name, seqDict[name].upper())
    #new = ['>%s\n%s\n' % (name, seqDict[name].upper().replace('-','N'))
            for name in orderlist]

    dn = cmn.lastName(fn).replace('.fa', '') + '_allUp.fa'
    cmn.write_file(''.join(new), dn)


