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


    #cmd = 'grep "Estimated Ln Prob of Data" %s/*/r*/*.log' % wdir
    cmd = 'ls %s/*/r*/structure*f' % wdir

    fns = cmn.cmd2lines(cmd)
    print(fns)

    outdir = 'harvest_%s' % wdir
    cmn.mkdir(outdir)

    cmd_dict = {}
    for fn in fns:
        #cov_3/structureK10/r0/structure.output_f
        K = cmn.find_between(fn, 'structureK', '/')
        rep = fn.split('/')[-2]
        dn = '%s/out_K%s_%s_f' % (outdir, K, rep)
        cmd = 'cp %s %s' % (fn, dn)
        try:
            cmd_dict[K].append(cmd)
        except KeyError:
            cmd_dict[K] = [cmd]

    for K in cmd_dict:
        cmds = cmd_dict[K]
        if len(cmds) < 3:
            print('insufficent replicates for K=%s, skip' % K)
            continue
        for cmd in cmds:
            cmn.run(cmd)
