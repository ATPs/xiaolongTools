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
        print("Usage: *.py coding.gff", file=sys.stderr)
        sys.exit()

    coding_indexes = []
    with open(fn) as fp:
        for line in fp:
            #scaffold1_cov14552_reverse	mitfi	trnM(cat)	23	89	1.809e-09	+	.
            items = line.strip().split()
            Range = list(map(int, items[3:5]))
            if items[6] == '-':
                j, i = Range
            else:
                i, j = Range
            indexes = list(range(i, j+1))
            coding_indexes += indexes


    dn = 'coding.indexes.pkl'
    cmn.pickle_write(set(coding_indexes), dn)



