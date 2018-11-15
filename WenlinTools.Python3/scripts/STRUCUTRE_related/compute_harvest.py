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
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        wdir = sys.argv[1]
    except:
        print("Usage: *.py wdir", file=sys.stderr)
        sys.exit()

    wdir = wdir.rstrip('/')
    cmd = 'make4harvest.py %s' % wdir
    cmn.run(cmd)

    cwd = os.getcwd()
    outdir = '%s/harvest_%s' % (cwd, wdir)

    #fake K1
    #cmd = 'cp /project/biophysics/Nick_lab/wli/sequencing/Phoebis/STRUCTURE/harvest_cov3_fewP/out_K1_r1_f %s' % outdir
    #cmn.run(cmd)

    rdir = 'harvest_report_%s' % wdir
    cmd = 'cd /home2/wli/local/structureHarvester\n'
    cmd += './structureHarvester.py --dir=%s --out=%s/%s --evanno' % (outdir, cwd, rdir)
    cmn.run(cmd)

    cmd = 'cd %s; cat evanno.txt' % rdir
    cmn.run(cmd)
