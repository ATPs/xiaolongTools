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


ddir = '/project/biophysics/Nick_lab/mtang/archive/step3_gatk'

wdir = sys.argv[1]

os.chdir(wdir)
vcf_list = [os.path.abspath(fn) for fn in cmn.cmd2lines('ls */*.vcf')]

spdirs = {'/'.join(each.split('/')[:-1]): each for each in vcf_list}

#deal with sample directory
kept_spdir_files = 'realigned_reads_step2.bam snp_step2.vcf$'.split()

for each in spdirs:
    print(each)
    wdir_label = cmn.lastName(each)
    dwdir = '%s/%s' % (ddir, wdir_label)
    if os.path.exists(dwdir):
        print('the destination directory has already exists! please check manually to choose which one to keep:')
        print('distination dir: %s' % dwdir)
        print('current dir: %s' % each) 
        print('\n')
        continue
    
    cmn.mkdir(dwdir)
    fbam = '%s/realigned_reads_step2.bam' % (each)
    if os.path.exists(fbam):
        cmd = 'cp %s/realigned_reads_step2.bam %s' % (each, dwdir)
    else:
        fbam = '%s/realigned_reads.bam' % each
        if os.path.exists(fbam):
            cmd = 'cp %s/realigned_reads.bam %s' % (each, dwdir)
        else:
            print('Error, can not find bam file!')
            
    cmn.run(cmd)
    #print cmd
    fvcf = spdirs[each]
    cmd = 'cp %s %s' % (fvcf, dwdir)
    cmn.run(cmd)
    #print cmd

    #try to clean

    cmd = 'rm -f %s; cd %s; ln -s %s/%s' % (fvcf, each, dwdir, cmn.lastName(fvcf))
    #cmn.run(cmd)

    cmd = 'rm -f %s/realigned_reads_step2.bam; cd %s; ln -s %s/realigned_reads_step2.bam' % (each, each, dwdir)
    #cmn.run(cmd)

