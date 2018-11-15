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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    adict = {}

    with open(fn) as fp:
        for line in fp:
            name, Nmapped, Ntotal = line.strip().split()
            sp = name.split('_')[0]
            Nmapped = int(Nmapped)
            Ntotal = int(Ntotal)
            try:
                current_mapped, current_total = adict[sp]
                adict[sp] = (current_mapped + Nmapped, current_total + Ntotal)
            except:
                adict[sp] = (Nmapped, Ntotal)

    print('\t'.join(['sp','Mapped','Total', 'Fraction']))
    for sp in adict:
        N1, N2 = adict[sp]
        print('%s\t%s\t%s\t%.2f%%' % (sp, N1, N2, (float(N1) / N2 * 100)))


