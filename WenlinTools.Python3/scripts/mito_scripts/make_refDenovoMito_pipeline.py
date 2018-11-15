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
        fqlist, fmitolist=sys.argv[1:]
    except:
        print("Usage: *.py fqlist refMitoList", file=sys.stderr)
        sys.exit()

    #ftemplate = '/work/biophysics/wli/Eudamine/wholeMito_run2/mito_denovo.template'
    ftemplate = '/project/biophysics/Nick_lab/wli/sequencing/scripts/mito_scripts/mito_refDenovo.template'
    template = cmn.txt_read(ftemplate)
    fqlist = cmn.file2lines(fqlist)
    groupDict = {}
    for fq in fqlist:
        fq = os.path.abspath(fq)
        ID = cmn.lastName(fq).split('_')[0]
        try:
            groupDict[ID].append(fq)
        except KeyError:
            groupDict[ID] = [fq]

    fmitolist = os.path.abspath(fmitolist)
    for sample in groupDict:
        fqlist = groupDict[sample]
        wdir = 'mitoRef_%s' % sample
        cmn.mkdir(wdir)
        os.chdir(wdir)
        cwd = os.getcwd()
        info = template.replace('[cwd]', cwd)
        cmn.write_lines(fqlist, 'fqlist')
        cmd = 'cat %s|xargs cat > ref_mito.fa; module add bwa; bwa index ref_mito.fa' % fmitolist
        cmn.run(cmd)

        bwa_commands = []
        if len(fqlist) == 3:
            R1 = [each for each in fqlist if '_R1' in each][0]
            R2 = [each for each in fqlist if '_R2' in each][0]
            single = [each for each in fqlist if '_single' in each][0]
            bwa_commands.append('bwa mem -M ref_mito.fa %s %s > paired.sam' % (R1, R2))
            bwa_commands.append('bwa mem -M ref_mito.fa %s > singleton.sam' % (single))
        elif len(fqlist) == 2:
            R1 = [each for each in fqlist if '_R1' in each][0]
            R2 = [each for each in fqlist if '_R2' in each][0]
            bwa_commands.append('bwa mem -M ../ref_mito.fa %s %s > paired.sam' % (R1, R2))
        else:
            print('don\'t know how to pair %s' % str(fqlist))

        info = info.replace('[bwa_commands]', '\n'.join(bwa_commands))

        dn = 'fD%s.job' % sample
        cmn.write_file(info, dn)

        os.chdir('..')




