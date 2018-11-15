import sys
import os


python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

#1. read in data
fns = cmn.getid(sys.argv[1])

falist = cmn.cmd2lines('ls *m2s.fa')

finished_maps = set([fn.replace('_m2s.fa', '.map') for fn in falist])

isGood = True

cmds = []
for fn in fns:
    label = cmn.lastName(fn)
    if label in finished_maps:
        continue

    isGood = False
    if 'MITO' in label:
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/map2fasta_mito.py %s' % fn
    else:
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/map2fasta.py %s' % fn
    cmds.append(cmd)

if isGood:
    print('Good news! everything looks good!')

else:
    cmds.append('')
    dn = 'm2fadd.cmds'
    cmn.write_lines(cmds, dn)

    print('Error!!!!!')
    print('There are still %s fasta missing' % (len(cmds) - 1))
    print('please use following command to submit jobs')
    print('\n>>> /work/biophysics/mtang/SNP_calling/scripts/submit_jobs.py %s [#node] m2fAdd -p 256GB\n' % dn)
    print('-p specifies the partition it submitted to')
    print('[#node] is the number of nodes and should be adjusted according to number of lines in %s' % dn)
    print('\n[IMPORTANT]Please run this check again upon the job completion.')

