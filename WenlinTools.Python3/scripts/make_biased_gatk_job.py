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
def group_by_species(fns):
    adict = {}
    for fn in fns:
        label = cmn.lastName(fn).split('_')[0]
        try:
            adict[label].append(fn)
        except:
            adict[label] = [fn]
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


if __name__=='__main__':
    #options=parse_options()

    fass = 'assembly_main.fa'
    ass_label = fass.replace('.fa', '')

    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/mapping/SNP_calling/2_gatk/template_gatk.job')
    template = template.replace('assembly_v0', ass_label)

    fns = cmn.cmd2lines('ls ../1_bwa_align/*.sam | grep -v _v0_')
    fns = [os.path.abspath(i) for i in fns]

    #good_set = set('LEP18259 3318 3303'.split())
    finished = cmn.cmd2lines("ls */*.vcf|cut -d '/' -f 2|cut -d '_' -f 1")

    #group spiecies
    group_dict = group_by_species(fns)

    for slabel in group_dict:
        if slabel in finished:
            print('skip finished ' + slabel)
            continue
        cmn.mkdir(slabel)
        os.chdir(slabel)

        fns = group_dict[slabel]
        f_sam = merge_sams(slabel, fns)

        cmd = template.replace('3377', slabel)
        #cmd = cmd.replace('--job-name=gatk', '--job-name=%s' % slabel)

        os.chdir('..')

        cmn.write_file(cmd, 'gatk%s.job' % slabel)
        cmn.run('sbatch gatk%s.job' % slabel)



