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
        fn = sys.argv[1]
        ref_scaf = sys.argv[2]
    except:
        print("Usage: *.py samfile scaf", file=sys.stderr)
        sys.exit()

    dn = 'filtered' + cmn.lastName(fn)
    dp = open(dn, 'w')
    with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                if ref_scaf in line:
                    dp.write(line)

            else:
                if line.strip().split()[2] == ref_scaf:
                    dp.write(line)

    dp.close()


