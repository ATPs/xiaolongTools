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


    rdict = {}
    start = False
    with open(fn) as fp:
        for line in fp:
            if 'Miss' in line:
                start = True
                continue
            elif 'Estimated Allele Frequencies in each cluster' in line:
                start = False

            if start:
                if line.strip() == '':
                    continue
                #63  15098E8    (0)   :  0.450 0.054 0.142 0.086 0.000 0.268
                items = line.strip().split()
                sp = items[1]

                rdict[sp] = '\t'.join(line.strip().split()[1:])
    new = []
    for sp in rdict:
        new.append(rdict[sp])
    print('\n'.join(new))

