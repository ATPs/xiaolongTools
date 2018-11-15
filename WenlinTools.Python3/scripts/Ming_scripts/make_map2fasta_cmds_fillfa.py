import sys
import os


python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

#1. read in data
fns = cmn.getid(sys.argv[1])
finished = [cmn.lastName(each) for each in cmn.getid(sys.argv[2])]

cmds = []
for fn in fns:
    label = cmn.lastName(fn)
    dn = label.replace('.map', '_m2s.fa')
    if dn in finished:
        print('skip finished %s' % dn)
        continue


    if 'MITO' in label:
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/map2fasta_mito.py %s' % fn
    else:
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/map2fasta.py %s' % fn
    cmds.append(cmd)


cmds.append('')
dn = 'm2f.cmds'
cmn.write_lines(cmds, dn)

print('After finish checking scaffold length, please use following command to submit jobs')
print('\n>>> /work/biophysics/mtang/SNP_calling/scripts/submit_jobs.py %s [#node] m2f -p 256GB\n' % dn)
print('-p specifies the partition it submitted to')
print('[#node] is the number of nodes and should be adjusted according to number of lines in %s' % dn)

