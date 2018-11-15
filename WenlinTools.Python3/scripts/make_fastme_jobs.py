import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tell_running_mode(fn):
    size = float(cmn.filesize(fn)) / 1024 / 1024 # M
    if size < 25:
        mode = 'RAxML'
    else:
        mode = 'ExaML'
    return mode
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wdir = os.path.abspath(sys.argv[1].rstrip('/'))

project = cmn.lastName(os.getcwd()).split('_')[0]

f_list = 'falist'

cmd = 'ls %s/* > %s' % (wdir, f_list)
cmn.run(cmd)

falist = [os.path.abspath(fn) for fn in cmn.file2lines(f_list)]

Njob = 3
fa_size = cmn.filesize(falist[0])/1024/1024

Njob = max(Njob, 50*fa_size/5000 + 1)

Ncores = 48 * Njob / 100
print('number of cores:', Ncores)
print('number of jobs:', Njob)

cmds = []
outdir = 'making_fastme_trees'
cmn.mkdir(outdir)
for fa in falist:
    cmd = 'cd %s; python /project/biophysics/Nick_lab/wli/sequencing/scripts/fasta2fastmeTree.py %s %s' % (outdir, fa, Ncores)
    cmds.append(cmd)

cmn.write_lines(cmds, 'fastme.cmds')

cmd = 'python /home2/wli/my_programs/submit_jobs.py fastme.cmds %s %s -p 256GB' % (Njob, project)
cmn.run(cmd)


