import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

jobs = [line.strip().split()[-1] for line in cmn.getid(sys.argv[1])]

fromDir = sys.argv[2].rstrip('/')#the dir ends with step3

cwd = os.getcwd()

cmn.mkdir('job_files')
cmn.mkdir('step3_gatk')

fromPdir = '/'.join(fromDir.split('/')[:-1])
cmn.run('ln -s %s/step2_bwa_mapping' % fromPdir)

fjobs = []
#1. copy the directory to current
for job in jobs:
    wdir = job[4:-4]
    current = '%s/%s' % (fromDir, wdir)
    cmd = 'cp -r %s step3_gatk' % current
    print('forking data for %s' % current)
    cmn.run(cmd)
    new = '%s/step3_gatk/%s' % (cwd, wdir)
    user = cmn.cmd2info('echo $USER').strip()
    user_label = user[0]

    fjob = '%s/job_files/%s' % (fromDir, job)
    info = cmn.txt_read(fjob)
    info = info.replace(fromDir, '%s/step3_gatk' % cwd)

    fjob = 'job_files/g%s%s.job' % (user_label, wdir)
    cmn.write_file(info, fjob)
    fjobs.append(cmn.lastName(fjob))

dn = 'forked_jobs.list'
cmn.write_lines(fjobs, dn)

