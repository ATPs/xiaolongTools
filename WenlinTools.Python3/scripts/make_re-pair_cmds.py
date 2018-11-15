#!/usr6/local/bin/python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #fns = cmn.cmd2lines('ls ../6.5_repair_ID/*.fq')
    fns = sys.argv[1:]

    pairDict = {}

    for fn in fns:
        key = '_'.join(cmn.lastName(fn).split('_')[:-1])
        #if '250' in key or '500' in key:
        #    print 'skip short lib: %s' % fn
        #    cmn.run('ln -s %s' % fn)
        #    continue

        try:
            pairDict[key].append(fn)
        except KeyError:
            pairDict[key] = [fn]

    cmn.mkdir('logs')
    cmds = []
    for key in pairDict:
        each = pairDict[key]
        if len(each) == 2:
            each.sort()
            cmd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/re-pair-reads_wenlin %s %s %s >& logs/%s_run.log &' % (each[0], each[1], key, key)
            cmds.append(cmd)
        else:
            print('cannot find pair for %s' % str(each))
            for iii in each:
                cmn.run('ln -s %s' % iii)

    cmds.append('wait')

    dn = 're-pair_auto.job'
    cmn.write_lines(cmds, dn)
