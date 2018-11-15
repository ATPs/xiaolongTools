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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    cmd = 'python /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/stepA6_throw_warning.py .'
    cmn.run(cmd)

    finished = set([each.split()[0].split('_')[0] for each in cmn.file2lines('checkSummary.report')])

    wdirs = cmn.cmd2lines('ls sampleRun_* -d ')

    cmds = []
    for wdir in wdirs:
        sp = wdir.split('_')[-1]
        if sp in finished:
            print('skip finished %s' % sp)
            continue
        cmds.append('cd %s; bash /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/error_detection.sh %s; cd ..;\n' % (wdir, sp))


    dn = 'Ecor_toRun.cmds'
    print('number of commands: %s' % len(cmds))
    print('commands is in %s' % dn)
    cmn.write_file(''.join(cmds), dn)

