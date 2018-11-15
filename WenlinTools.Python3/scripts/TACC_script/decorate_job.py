#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
#main

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import sys
python_lib='/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

if __name__=='__main__':
    import cmn,os,sys
    k=len(sys.argv)

    if k == 1:
        print('usage:make_job.py fn [-n 4 -p 128G -t 24]')
        sys.exit()

    fn=sys.argv[1]
    key_cmd = cmn.txt_read(fn)

    node = '1'
    part = 'normal'
    time_hour = '48'
    for i, arg in enumerate(sys.argv):
        if arg == '-n':
            node = sys.argv[i+1]
        elif arg == '-p':
            part = sys.argv[i+1]
        elif arg == '-t':
            time_hour = sys.argv[i+1]


    cwd=os.getcwd()

    aa=cmn.txt_read('/work/00412/mtang/sequencing/scripts/slurm.job')

    aa = aa.replace('NODE', node)
    aa = aa.replace('PART', part)
    aa = aa.replace('TIME_HOUR', time_hour)

    #aa+="#$ -pe %sway %s\n\n\n"  % (cpu, cpu)

    aa+='cd %s\n\n' % cwd

    aa+=key_cmd+'\n'

    print(aa)

