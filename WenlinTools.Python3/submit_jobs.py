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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


import os
import cmn

def submit_job(dn, cores, pre, i, partition):
    tmpdir, dn2 = dn.split('/')
    cmd = 'source /home2/wli/.bash_profile;'
    incmd = 'master_run_jobs.py %s %s' % (dn, cores)
    cmd += 'make_job.py "%s" -p %s > %s/W%s_%s.job;' % (incmd, partition, tmpdir, pre, i)
    cmd += 'cd %s;' % tmpdir
    cmd += 'sbatch W%s_%s.job' % (pre, i)
    os.system(cmd)


if __name__ == '__main__':
    try:
        cmds = cmn.file2lines(sys.argv[1])
        N = int(sys.argv[2])
    except:
        print('usage: *.py cmd_file split_parts [prefix4spilt,def:x] [-p super] [-t 24]', file=sys.stderr)
        sys.exit()

    try:
        pre = sys.argv[3]
    except:
        pre = 'x'

    cores = 32 # in biohpc, it is 32

    partition = 'super'
    for i, arg in enumerate(sys.argv):
        if arg == '-p':
            partition = sys.argv[i+1]

    tmpdir = '%s_files' % pre
    cmn.mkdir(tmpdir)
    bundle = len(cmds) / N

    for i in range(N):
        cmd = cmds[i*bundle: (i+1)*bundle]
        if i == N - 1:
            cmd += cmds[(i+1)*bundle:]
        dn = '%s/%s_%s' % (tmpdir, pre, i)
        cmn.write_file('\n'.join(cmd), dn)
        submit_job(dn, cores, pre, i, partition)

    cwd = os.getcwd()
    d_stat = '%s/stat.info' % tmpdir
    info = []
    info.append('commands are from: %s/%s' % (cwd,sys.argv[1]))
    info.append('split into %s jobs' % sys.argv[2])
    info.append('submitted to queue: %s' % partition)
    cmn.write_file('\n'.join(info), d_stat)
