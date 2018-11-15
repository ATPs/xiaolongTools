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
import numpy as np


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fns=sys.argv[1:]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    lengths = []
    for fn in fns:
        with open(fn) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 1:
                    lengths.append(len(line.strip()))


    print(fn, np.mean(lengths), np.std(lengths, ddof=1))




