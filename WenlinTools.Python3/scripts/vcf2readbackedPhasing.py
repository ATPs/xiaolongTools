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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fvcf, fbam = list(map(os.path.abspath, sys.argv[1:]))
    except:
        print("Usage: *.py *.vcf *.bam", file=sys.stderr)
        print("Generate the command to phase vcf", file=sys.stderr)
        sys.exit()


    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/template_readbackedPhasing.cmds')
    dnlabel = cmn.lastName(fvcf).replace('.vcf', '')
    outdir = '%s_wdir' % dnlabel
    cmn.mkdir(outdir)
    os.chdir(outdir)

    cwd = os.getcwd()

    cmd = 'ln -s %s' % fvcf
    cmn.run(cmd)
    fvcf = cmn.lastName(fvcf)

    cmd = 'ln -s %s' % fbam
    cmn.run(cmd)
    fbam = cmn.lastName(fbam)

    info = template.replace('[input.bam]', fbam)

    info = info.replace('[input.vcf]', fvcf)

    info = info.replace('[preprocessing_cmds]', 'cd %s' % cwd)


    dn = 'phased_%s.vcf' % cmn.lastName(fvcf)
    info = info.replace('[output.vcf]', dn)

    info = info.replace('[windowSize]', '2000')

    flog = '%s.log' % dnlabel
    info = info.replace('[log_file]', flog)

    dn = 'ph%s.job' % (dnlabel)
    cmn.write_file(info, dn)

    print('bash %s/%s' % (cwd, dn))
