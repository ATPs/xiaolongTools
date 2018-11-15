import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os
import time

def get_current_jobs(label, user):
    cmd = 'squeue| grep %s| grep g%s|wc -l' % (user, label)
    N = cmn.cmd2info(cmd).split()[0]
    N = int(N)
    return N


fn = 'forked_jobs.list'

jobs = cmn.getid(fn)

cores = int(sys.argv[1])

user = cmn.cmd2info('echo $USER').strip()
user_label = user[0]

currentN = get_current_jobs(user_label, user)

os.chdir('job_files')

todo = list(jobs)

while(len(todo) != 0):
    fjob = todo[0]
    currentN = get_current_jobs(user_label, user)
    print(currentN)
    if currentN < cores:
        #submit
        cmd = 'sbatch %s' % fjob
        cmn.run(cmd)
        todo.remove(fjob)
    else:
        time.sleep(300)
