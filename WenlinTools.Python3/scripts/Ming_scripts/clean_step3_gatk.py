import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os


#files to keep
#1. job_files (but remove every job*err and job*out)
#2. in each sample directory, keep 
#   assembly_selfref_v2.fa, realigned_reads_step2.bam, 5660_snp_step2.vcf, 


wdir = sys.argv[1]

os.chdir(wdir)
print('cd %s' % wdir)

fns = cmn.cmd2lines('ls ;')

#keep some files
kept = 'job_files'.split()
for each in kept:
    try:
        fns.remove(each)
    except:
        continue

#clean job_files
cmn.run('rm job_files/job* 2> /dev/null')

#deal with sample directory
spdirs = set([line.strip().split('/')[0] for line in cmn.cmd2lines('ls */*.vcf')])
kept_spdir_files = 'assembly_selfref_v2.fa$ realigned_reads_step2.bam snp_step2.vcf$'.split()

for each in spdirs:
    try:
        fns.remove(each)
    except:
        continue
    
    #clean each directory
    cmd = 'ls %s/*|' % each
    for spfn in kept_spdir_files:
        cmd += 'grep -v %s|' % spfn
    cmd += 'xargs rm -r' 
    print(cmd)
    #cmn.run(cmd)

for each in fns:
    cmd = 'rm -r %s' % each
    #cmn.run(cmd)
    print(cmd)
