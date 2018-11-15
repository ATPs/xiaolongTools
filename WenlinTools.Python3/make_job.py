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
python_lib='/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)



if __name__=='__main__':
    import cmn,os,sys
    k=len(sys.argv)

    if k == 1:
        print('usage:make_job.py "key_command" [-n 4 -p 128G]')
        sys.exit()

    key_cmd=sys.argv[1]

    node = '1'
    part = '256GB'
    time_hour = '200'
    memLimit = None
    for i, arg in enumerate(sys.argv):
        if arg == '-n':
            node = sys.argv[i+1]
        elif arg == '-p':
            part = sys.argv[i+1]
        elif arg == '-t':
            time_hour = sys.argv[i+1]
        elif arg == '-m':
            memLimit = sys.argv[i+1]


    adding = []
    if memLimit != None:
        adding += '#SBATCH --mem %s\n' % memLimit

    cwd=os.getcwd()

    aa=cmn.txt_read('/home2/wli/template/slurm.job')

    aa = aa.replace('NODE', node)
    aa = aa.replace('PART', part)
    aa = aa.replace('TIME_HOUR', time_hour)
    aa = aa.replace('[variables]', ''.join(adding))

    #aa+="#$ -pe %sway %s\n\n\n"  % (cpu, cpu)

    aa+='cd %s\n\n' % cwd

    aa+=key_cmd+'\n'

    print(aa)

