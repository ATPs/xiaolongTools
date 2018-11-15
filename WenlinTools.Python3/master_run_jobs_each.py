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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


from multiprocessing import Pool, cpu_count
import os
import cmn

#your function to run your job
def run_job(cmd):
    #print 'running command:%s' % cmd
    os.system(cmd)
    return 'finish running: %s' % cmd


if __name__ == '__main__':
    try:
        cmds = cmn.file2lines(sys.argv[1])
        cores = int(sys.argv[2])
    except:
        print('usage: *.py cmd_file split_parts', file=sys.stderr)
        sys.exit()

    #detect cores by the machine
    #cores = cpu_count()
    #print 'running jobs using %s cpus' % cores
    pool = Pool(processes=cores)              # start 4 worker processes

    result_list = []
    for line in cmds:
        cmd = line.strip()
        process = pool.apply_async(run_job, [cmd])
        result_list.append(process)

    #pull out the result
    for process in result_list:
        #print process.get()
        process.get()
