#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.clw", file=sys.stderr)
        sys.exit()

    new = []
    with open(fn) as fp:
        for line in fp:
            exon, sp, seq = line.strip().split()
            sp = sp.split('.')[0]
            new.append('>%s_%s\n%s\n' % (sp, exon, seq))

    dn = cmn.lastName(fn).replace('.sum', '') + '.fa'
    cmn.write_file(''.join(new), dn)


