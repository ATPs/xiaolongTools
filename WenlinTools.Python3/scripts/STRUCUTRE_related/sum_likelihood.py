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
        wdir = sys.argv[1]
    except:
        print("Usage: *.py t100_highQ", file=sys.stderr)
        sys.exit()


    cmd = 'grep "Estimated Ln Prob of Data" %s/*/r*/*.log' % wdir

    lines = cmn.cmd2lines(cmd)
    print('\n'.join(lines))

    rdict = {}
    countK = {}
    for line in lines:
        K = cmn.find_between(line, 'structureK', '/')
        K = int(K)
        lnL = float(line.strip().split()[-1])
        try:
            rdict[K].append(lnL)
        except KeyError:
            rdict[K] = [lnL]

        try:
            countK[K] += 1
        except:
            countK[K] = 1

    keys = list(rdict.keys())
    keys.sort()
    for K in keys:
        rep = countK[K]
        if rep == 1:
            print('K=%s, rep=%s, lnL=%s ' % (K, countK[K], np.mean(rdict[K])))
        else:
            print('K=%s, rep=%s, lnL=%s (STD:%s, max:%s)' % (K, countK[K], np.mean(rdict[K]), np.std(rdict[K], ddof=1), max(rdict[K])))




