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
        fn=sys.argv[1]
    except:
        print("Usage: *.py gap_see", file=sys.stderr)
        sys.exit()

    pop_map = {}
    fpops = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/general_info/P*IDs')
    for fpop in fpops:
        alist = cmn.getid(fpop)
        popname = cmn.lastName(fpop)[1:-3]
        for sp in alist:
            pop_map[sp] = popname

    #15101E04_snp.codeVcf_father 13629191 1986713 0.145768960168
    adict = {}
    for line in cmn.file2lines(fn):
        items = line.split()
        sp = items[0].split('_')[0]
        gapF = float(items[-1])
        try:
            pop = pop_map[sp]
        except:
            continue

        try:
            adict[pop].append(gapF)
        except:
            adict[pop] = [gapF]

    for pop in adict:
        print(pop, '%.2f' % np.mean(adict[pop]))

    for pop in adict:
        print(pop, max(adict[pop]))



