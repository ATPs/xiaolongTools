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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import numpy as np

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()

    gap_chars = ['N', '-']
    alist = []
    totalN = 0
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line.strip()[1:]
            else:
                seq = line.strip()
                totalN += len(seq)
                Ngap = 0
                for char in gap_chars:
                    Ngap += seq.count(char)
                value = float(Ngap)/len(seq)
                alist.append(value)

    print(fn, totalN, np.mean(alist))




