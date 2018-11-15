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
        Ncores = int(sys.argv[2])
    except:
        print("Usage: *.py fa Ncores", file=sys.stderr)
        sys.exit()

    fphylip = cmn.lastName(fn) + '.phylip'
    cmd = 'python /project/biophysics/Nick_lab/wli/sequencing/scripts/fasta2phylip4fastme.py %s' % fn
    cmn.run(cmd)

    dn = fphylip + '.fastme.dist'
    cmd = '/home2/wli/local/bin/fastme -i %s -O %s -d1 -T %s' % (fphylip, dn, Ncores)
    cmn.run(cmd)




