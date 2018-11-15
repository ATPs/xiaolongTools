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
        indir, label = sys.argv[1:3]
    except:
        print("Usage: *.py ../1_process_alignment/noGap_splits noGap", file=sys.stderr)
        sys.exit()

    fns = [os.path.abspath(i) for i in cmn.cmd2lines('ls %s/*' % indir)]

    wdir = 'split_run_%s' % label
    cmn.mkdir(wdir)
    os.chdir(wdir)
    cmn.mkdir('job_files')
    for count, fn in enumerate(fns):
        cmn.run('ln -s %s' % fn)
        fn_new = cmn.lastName(fn)
        cmd = 'rm *%s_%s; /home2/wli/local/RAxML/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -s %s -n %s_%s -p 7112 -T 48 ' % (label, count, fn_new, label, count)
        dn = 'job_files/sg%s.job' % count
        cmn.run('/home2/wli/my_programs/make_job.py "%s" -p 256GB -t 33 > %s' % (cmd, dn))




