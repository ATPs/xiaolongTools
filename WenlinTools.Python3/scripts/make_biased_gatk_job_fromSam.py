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
    try:
        fsam, fass = list(map(os.path.abspath, sys.argv[1:]))
    except:
        print('usage: *.py fsam fass', file=sys.stderr)
        sys.exit()

    cmd = 'module add samtools; samtools faidx %s' % fass
    cmn.run(cmd)
    cmd = 'module add picard/1.117; java -jar $PICARD/CreateSequenceDictionary.jar R=%s O=%s.dict' % (fass, fass[:-3])
    cmn.run(cmd)

    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/template_gatk_bias_fromSam.job')
    template = template.replace('[WL_ref]', fass)
    template = template.replace('[INPUT.sam]', fsam)

    sampleId = cmn.lastName(fsam).replace('highQ_', '').split('_')[0]

    dnlabel = '%s_%s' % (cmn.lastName(fsam).replace('.sam', ''), cmn.lastName(fass).replace('.fa', ''))
    cmn.mkdir(dnlabel)
    os.chdir(dnlabel)

    cwd = os.getcwd()
    pre_cmds = 'cd %s\n' % cwd
    template = template.replace('5642', sampleId)
    template = template.replace('[WL_preprossing]', pre_cmds)

    cmn.write_file(template, 'gatk%s.job' % sampleId)
    #cmn.run('sbatch gatk%s.job' % sampleId)



