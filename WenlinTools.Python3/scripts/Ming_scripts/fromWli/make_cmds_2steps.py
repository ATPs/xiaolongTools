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

import cmn
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def separate_by_label(fns):
    adict = {}
    for fn in fns:
        label = '_'.join(cmn.lastName(fn).split('_')[:-1])
        try:
            adict[label].append(fn)
        except:
            adict[label] = [fn]
    return adict


def separate_by_pair(label, fns):
    #paired = [i for i in fns if ('_paired' in i) ]
    paired = [i for i in fns if ('_paired' in i) or ('_R' in i) or ('_pair1' in i) or ('_pair2' in i)]
    paired.sort()
    if len(paired) != 2:
        print('error: wrong number of pairs as %s' % str(paired))
        print('from: %s' % str(fns))
        print('need to change the label criterion')
        sys.exit()

    unpaired = set(fns) - set(paired)

    #parse each files
    newPaired = []
    for fn in paired:
        cmn.run('unlink %s; ln -s %s' % (cmn.lastName(fn), fn))
        newPaired.append(cmn.lastName(fn))

    singleFn = '%s_single.fq' % label
    if cmn.filexist:
        cmn.run('rm %s' % singleFn)
    for fn in unpaired:
        cmn.run('cat %s >> %s' % (fn, singleFn))

    return newPaired, singleFn


def make_bwa_cmds(fqs, f_ass):
    group_dict = separate_by_label(fqs)

    fsams = []
    for plabel in group_dict:
        print('processing lib %s' % plabel)
        each = group_dict[plabel]
        #also parse the files inside this function
        #return the file name after parsing
        paired, unpaired = separate_by_pair(plabel, each)
        #label = '%s_%s' % (plabel, ass_label)
        index_label = cmn.lastName(f_ass).replace('.fa', '')
        cmd = ''
        cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 48 -M %s %s %s > %s_paired.sam;\n' % (index_label, paired[0], paired[1], plabel)
        cmd += '/home2/wli/local/bwa-0.7.12/bwa mem -t 48 -M %s %s > %s_single.sam;\n' % (index_label, unpaired, plabel)
        fsams.append('%s_paired.sam' % plabel)
        fsams.append('%s_single.sam' % plabel)

    return cmd, fsams



def merge_sams(label, fns):
    dn = '%s.sam' % label

    print('merging files: %s into %s' % (str(fns), dn))

    if cmn.filexist(dn):
        cmn.run('rm ' + dn)

    cmn.run('cp %s %s' % (fns[0], dn))

    fp_dn = open(dn,"a")
    for fn in fns[1:]:
        fp = open(fn)
        for line in fp:
            if line[0] != "@" and line[0] != "[" and line.split()[2] != "*":
            #if line[0] != "@":
                fp_dn.write(line)
        fp.close()

    fp_dn.close()
    return dn



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py filelist", file=sys.stderr)
        sys.exit()

    #step 1 including 1. making the pileup file, 2. correct bias, 3. sam map again 4. snp call till last step
    #step 2: just the snp call using multiple CPUs

    #only put 15 jobs in a node to avoid memmory problem
    f_ass = '/project/biophysics/Nick_lab/wli/sequencing/Nick_request/Heli_map/SNP_calling/2_gatk/assembly_v2.fa'
    #should not run for reference genome; this should just build by original snp call
    ref_sp = '3935'

    step1dir = 'step1_cmds'
    cmn.mkdir(step1dir)

    #fns = cmn.getid(fn)
    fns = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Nick_request/Heli_map/SNP_calling/2_gatk/*/realigned_reads.bam')
    finished_pileups = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Nick_request/Heli_map/SNP_calling/2_gatk/*/*.pileup')

    finished = set([cmn.lastName(i).split('_')[0] for i in cmn.cmd2lines('ls */*.vcf')])
    step1_finished = set([cmn.lastName(i).split('/')[0] for i in cmn.cmd2lines('ls */realigned_reads.bam')])

    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/myAnalysis/clean_ref_bias/1_gatk_runs/step1_job.template')
    template2 = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/myAnalysis/clean_ref_bias/1_gatk_runs/step2_job.template')

    cmn.mkdir('job_files')

    step1_jobs = []
    step2_jobs = []

    good_sps = '5658 5950 15105C07 5449 5962'.split()

    for fn in fns:
        items = fn.split('/')
        sp = items[-2]
        if sp == ref_sp:
            continue

        if sp in finished:
            print('skip finished %s' % sp)
            continue

        #if sp not in good_sps:
        #    continue
        #cmn.run('rm -r %s' % sp)
        cmn.mkdir(sp)

        os.chdir(sp)

        step1cmds = []

        #step1.1, make the pileup
        f_pileup = fn.replace('/realigned_reads.bam', '/realigned.pileup')
        if f_pileup not in finished_pileups:
            cmd = '/home2/wli/local/samtools-1.2/samtools mpileup -A -B -x -d 999999 %s -f %s > %s' % (fn, f_ass, f_pileup)
            step1cmds.append(cmd)

        #step1.2, make the unbias ref
        #f_unbias = fn.replace('realigned_reads.bam', 'assembly_selfref_v2.fa')
        #if not cmn.filexist(f_unbias):
        cmd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/change2dominateForm_withRef.py %s %s' % (f_pileup, f_ass)
        step1cmds.append(cmd)

        #after getting the ambiguos genome, use it as fn in the following
        #fn = f_unbias
        fn = 'assembly_selfref_v2.fa'
        asslabel = fn.strip().split()[-1].split('.')[0]

        #cmn.run('ln -s %s' % fn)
        #fn = cmn.lastName(fn)

        fqs = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Nick_request/Heli_map/SNP_calling/1_bwa_align/%s*.fq' % sp)
        fqs += cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Nick_request/Heli_map/SNP_calling/1_bwa_align/%s*.fastq' % sp)
        if len(fqs) == 0:
            print('can not find fastq for %s!' % sp)
            os.chdir('..')
            continue

        cmd, fsams = make_bwa_cmds(fqs, fn)
        dn = 'sam_filelist'
        cmn.write_lines(fsams, dn)

        info = template.replace('assembly_selfref', asslabel)
        info = info.replace('[bwa_cmds]', cmd)
        info = info.replace('[WL_sam_filelist]', dn)
        info = info.replace('[WL_preprocessing]', '\n'.join(step1cmds))

        #make snp call cmds
        #f_sam = merge_sams(sp, fsams)

        info = info.replace('5328', sp)
        info = info.replace('[WL_cwd]', os.getcwd())

        info2 = template2.replace('assembly_selfref', asslabel)
        info2 = info2.replace('5328', sp)
        info2 = info2.replace('[WL_cwd]', os.getcwd())

        os.chdir('..')
        fjob = 'job_files/s1_%s.job' % sp
        cmn.write_file(info, fjob)
        cmn.run('cd job_files; sbatch s1_%s.job' % sp)

        if sp not in step1_finished:
            step1_jobs.append(fjob)

        fjob = 'job_files/s2_%s.job' % sp
        cmn.write_file(info2, fjob)
        step2_jobs.append(fjob)

        #cmn.run('cd job_files; sbatch sg%s.job' % sp)

    info = ['bash %s\n' % each for each in step1_jobs]
    cmn.write_file(''.join(info), 'step1todo.cmds')

    info = ['sbatch %s\n' % each for each in step2_jobs]
    cmn.write_file(''.join(info), 'step2todo.cmds')


