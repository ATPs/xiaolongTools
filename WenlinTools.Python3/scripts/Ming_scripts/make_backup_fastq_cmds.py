import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

cwd = os.getcwd()
wdir = os.path.abspath(sys.argv[1].rstrip('/'))

os.chdir(wdir)

fastqs = cmn.cmd2lines('ls *q')

groups = {}

for fq in fastqs:
    sp = fq.split('_')[0]
    try:
        groups[sp].append(fq)
    except:
        groups[sp] = [fq]

newlist = []
oldlist = []
for sp in groups:
    fqs = groups[sp]
    if len(fqs) == 2:
        #no single
        fpaired = fqs
        fpaired.sort()

        cmd = 'python /project/biophysics/Nick_lab/mtang/archive/update_ready_fastq_paired.py %s' % ' '.join(fpaired)
        newlist.append('%s_R1.fastq' % sp)
        newlist.append('%s_R2.fastq' % sp)

        print(cmd)
        oldlist += fqs
    elif len(fqs) == 3:    
        fsingle = ['%s/%s' % (wdir, fn) for fn in fqs
                if 'single' in fn][0]
    
        fpaired = ['%s/%s' % (wdir, fn) for fn in fqs
                if 'single' not in fn]
        fpaired.sort()

        cmd = 'python /project/biophysics/Nick_lab/mtang/archive/update_ready_fastq.py %s singleton' % fsingle
        print(cmd)
        newlist.append('%s_singleton.fastq' % sp)

        cmd = 'python /project/biophysics/Nick_lab/mtang/archive/update_ready_fastq_paired.py %s' % ' '.join(fpaired)
        newlist.append('%s_R1.fastq' % sp)
        newlist.append('%s_R2.fastq' % sp)

        print(cmd)
        oldlist += fqs
    elif len(fqs) == 4:
        fsingle = ['%s/%s' % (wdir, fn) for fn in fqs
                if '_unpaired' in fn]
    
        fpaired = ['%s/%s' % (wdir, fn) for fn in fqs
                if '_paired' in fn]
        fpaired.sort()

        for fn in fsingle:
            cmd = 'python /project/biophysics/Nick_lab/mtang/archive/update_ready_fastq.py %s singleton' % fn
            print(cmd)
        newlist.append('%s_singleton.fastq' % sp)

        cmd = 'python /project/biophysics/Nick_lab/mtang/archive/update_ready_fastq_paired.py %s' % ' '.join(fpaired)
        newlist.append('%s_R1.fastq' % sp)
        newlist.append('%s_R2.fastq' % sp)
        print(cmd)
        oldlist += fqs

    
    else:
        print('Error! wrong number for %s! expect 3 libs! skip' % sp)
        continue

cmd = 'cd %s' % wdir
print(cmd)
ddir = '/project/biophysics/Nick_lab/mtang/archive/fastq_libs'

for fn in oldlist:
    cmd = 'rm -f %s' % fn
    print(cmd)

for fn in newlist:
    cmd = 'ln -s %s/%s' % (ddir, fn)
    print(cmd)

