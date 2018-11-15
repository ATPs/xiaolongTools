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
        fn=sys.argv[1]
    except:
        print("Usage: *.py fqlist", file=sys.stderr)
        sys.exit()

    #ftemplate = '/work/biophysics/wli/Eudamine/wholeMito_run2/mito_denovo.template'
    ftemplate = '/project/biophysics/Nick_lab/wli/sequencing/scripts/mito_scripts/mito_denovo.template'
    template = cmn.txt_read(ftemplate)
    fqlist = cmn.file2lines(fn)
    groupDict = {}
    for fq in fqlist:
        ID = cmn.lastName(fq).split('_')[0]
        try:
            groupDict[ID].append(fq)
        except KeyError:
            groupDict[ID] = [fq]


    for sample in groupDict:
        fqlist = groupDict[sample]
        #fqlist = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/tmp_link_fastq/%s*.fastq' % sample)
        #fqlist = cmn.cmd2lines('ls /work/biophysics/wli/workspace/filtered_6313*q')
        wdir = 'mitoD_%s' % sample
        cmn.mkdir(wdir)
        os.chdir(wdir)
        cwd = os.getcwd()
        info = template.replace('[cwd]', cwd)
        info = info.replace('[fq_files]', ' '.join(fqlist))
        info = info.replace('[sample]', sample)

        #prepare quake infiles
        fqlist_local = []
        for fq in fqlist:
            cmn.run('ln -s ' + fq)
            fqlist_local.append(cmn.lastName(fq))
        cmn.write_lines(fqlist_local, 'fqlist')
        cmn.run('ln -s fqlist infiles')

        #make fq2fa comand
        quake_fqlist = [each.replace('.fastq', '.cor.fastq') for each in fqlist_local]
        fq2fa_cmds = ['rm %s.fa 2> /dev/null' % sample]
        for fq in quake_fqlist:
            cmd = 'fq2fa %s >> %s.fa;' % (fq, sample)
            fq2fa_cmds.append(cmd)

        cmn.write_lines(quake_fqlist, 'fqlist.cor')
        cmd = '\n'.join(fq2fa_cmds)
        info = info.replace('[fq2fa_commands]', cmd)

        noWolba_fqlist = ['noWolb_%s' % each for each in quake_fqlist]
        info = info.replace('[noWolba_fastq]', ' '.join(noWolba_fqlist))
        R1 = [each for each in noWolba_fqlist if '_R1' in each]
        R2 = [each for each in noWolba_fqlist if '_R2' in each]
        noWolba_pair = R1 + R2
        info = info.replace('[noWolba_pair]', ' '.join(noWolba_pair))

        dn = 'mD%s.job' % sample
        cmn.write_file(info, dn)

        cmd = 'sbatch %s' % dn
        cmn.run(cmd)

        os.chdir('..')




