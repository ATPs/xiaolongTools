import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

search_dirs = ['/project/biophysics/Nick_lab/mtang/archive/fastq_libs',
                '/project/biophysics/Nick_lab/mtang/preprocessing/*']

def search_for_old_fastq(label, wdirs):
    global selfDir
    alist = []
    for wdir in wdirs:
        cmd = 'ls %s/%s*q 2> /dev/null' % (wdir, label)
        alist += [line for line in cmn.cmd2lines(cmd)
                if selfDir not in line]
    return alist        

def remove_duplication(alist):
    stat_dict = {}
    dup = []
    for fn in alist:
        cmd = 'wc -l %s' % fn
        N = int(cmn.cmd2info(cmd).strip().split()[0])
        if N in stat_dict:
            dup.append(fn)
        else:
            stat_dict[N] = fn
    return list(stat_dict.values()), dup            
        
def group_fastq(fns):
    adict = {}
    for fn in fns:
        key = cmn.lastName(fn).split('_')[0]
        try:
            adict[key].append(fn)
        except KeyError:
            adict[key] = [fn]
    return adict


#~~~~~~~~main~~~~~~~~~~~~~~#    
wdir = sys.argv[1].rstrip('/')
selfDir = os.path.abspath(wdir)

fastqs = cmn.cmd2lines('ls %s/*q' % wdir)
print(fastqs)

log_info = []
outdir = '%s/fastq_thisBatch' % wdir
cmn.mkdir(outdir)

hasCombined = False
for fastq in fastqs:
    label = '.'.join(cmn.lastName(fastq).split('.')[:-1])
    print('processing %s' % label)
    old_fastqs = search_for_old_fastq(label, search_dirs)
    dn = '%s/%s' % (outdir, cmn.lastName(fastq))
    if not os.path.exists(dn):
        cmn.run('mv %s %s' % (fastq, dn))
    else:
        #has processed this lib
        continue

    hasCombined = True
    if len(old_fastqs) == 0: # no old data
        print('no old libs found for %s' % label)
        cmn.run('ln -s %s' % dn)
    else:#has old data
        print('combining old libs for %s' % label)
        old_fastqs, dup_fastqs = remove_duplication(old_fastqs)
        cmn.run('cp %s %s' % (dn, wdir))
        log_info.append('%s\t%s\n' % (label, dn))
        comb_fn = '%s/%s' % (wdir, cmn.lastName(dn))
        for old_fastq in old_fastqs:
            cmn.run('cat %s >> %s' % (old_fastq, comb_fn))
            log_info.append('%s\t%s\n' % (label, old_fastq))

if hasCombined:    
    cmn.write_file(''.join(log_info), '%s/combined_libs.log' % wdir)    

#make statistics for data amount

fastq_groups = group_fastq(fastqs)

new = []
for key in fastq_groups:
    fns = fastq_groups[key]
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/check_fastq_size.py %s %s' % (key, ','.join(fns))
    new.append(cmd)

new.append('')
cmn.write_lines(new, 'fastq_amount.cmds')    

