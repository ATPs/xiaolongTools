#!/usr/bin/env python

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
    #fns = cmn.cmd2lines('ls ../0_libs/*/*.fq')
    #fns += cmn.cmd2lines('ls ../0_libs/*/*.fastq')
    #fns = cmn.file2lines('../fqlist')
    #fns = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/tmp_link_fastq/*q')
    fns = cmn.file2lines(sys.argv[1])

    #skip_list = set(['5316', '5721'])

    gdict = {}
    for fn in fns:
        #items = fn.split('/')
        #sp = items[-2]
        sp = cmn.lastName(fn).split('_')[0]
        try:
            gdict[sp].append(fn)
        except:
            gdict[sp] = [fn]

    formatcmds = '\n\n'
    for sp in gdict:
        #if sp in skip_list:
        #    continue

        cmd = ''
        fns = gdict[sp]
        for fn in fns:
            cmd += 'fq2fa %s >> %s.fa; ' % (fn, sp)


        formatcmds += 'formatdb -p F -i %s.fa & \n' % sp
        print(cmd)

    #print formatcmds
    #print '\nwait\n'
