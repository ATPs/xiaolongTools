#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    alist = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
        alist.append(defline)
    return adict, alist



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py assembly", file=sys.stderr)
        sys.exit()


    seqDict, order_scafs = read_fa(fn)

    #write in the same dict as assembly
    olabel = cmn.lastName(fn)
    dnlabel = '.'.join(olabel.split('.')[:-1]) + '_scaf.header'
    dn = fn.replace(olabel, dnlabel)

    fdn = open(dn, 'w')

    for name in order_scafs:
        shortname = name.split()[0]
        length = len(seqDict[name])
        for i in range(length):
            line = '%s\t%s\n' % (shortname, (i+1))
            fdn.write(line)

    fdn.close()
