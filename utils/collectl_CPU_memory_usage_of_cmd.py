#!/usr/bin/env python
# encoding: utf-8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import os, sys, re
import logging
import argparse
import collections
import datetime
import time
import threading
import uuid
import signal
import subprocess
import psutil
import pandas as pd

logging.basicConfig(stream=sys.stderr, level=logging.INFO)
logger = logging.getLogger(__name__)

class Worker(threading.Thread):
    def __init__(self, shell):
        self.shell = shell
        self.proc = None
        threading.Thread.__init__(self)
    def start(self):
        self.proc = subprocess.Popen(self.shell, shell=True)
    def stop(self):
        time.sleep(3)
        parent = psutil.Process(self.proc.pid)
        children = parent.children(recursive=True)
        for p in children:
            if p.is_running():
                os.kill(p.pid, signal.SIGKILL)
        if parent.is_running():
            os.kill(parent.pid, signal.SIGKILL)
        



"""
Formatting from collectl -sZ

# PROCESS SUMMARY (counters are /sec)
# PID  User     PR  PPID THRD S   VSZ   RSS CP  SysT  UsrT Pct  AccuTime  RKB  WKB MajF MinF Command
    1  root     20     0    0 S   32M  552K  2  0.00  0.00   0  00:04.20    0    0    0    0 /sbin/init
    2  root     20     0    0 S     0     0  3  0.00  0.00   0  00:02.56    0    0    0    0 kthreadd
    3  root     RT     2    0 S     0     0  0  0.00  0.00   0  41:48.40    0    0    0    0 migration/0
    4  root     20     2    0 S     0     0  0  0.00  0.00   0  00:36.52    0    0    0    0 ksoftirqd/0

"""
# collectl = 'collectl'
# interval = 1
# cmd_job = 'python /home1/xc278/w/GitHub/perGeno/src/perGeno.py -g /gpfs/gpfs/project1/jx76-001/xc/20200812perGeno/20201216Test/GRCh38.p13.genome.fa.gz -p /gpfs/gpfs/project1/jx76-001/xc/20200812perGeno/20201216Test/gencode.v36.pc_translations.fa.gz -f /gpfs/gpfs/project1/jx76-001/xc/20200812perGeno/20201216Test/gencode.v36.chr_patch_hapl_scaff.annotation.gtf.gz -m /gpfs/gpfs/project1/jx76-001/xc/20200812perGeno/AFmostCommon/eas_adj.csv.gz -t 12 -o  /gpfs/gpfs/project1/jx76-001/xc/20200812perGeno/20201216Test/GENCODE12'
# PROGS_OF_INTEREST = ['perGeno']


def compute_GB(memory_val):

    r = re.search("^(\d+)([MKG])$", memory_val)

    if r:
        numG = int(r.group(1))
        metric = r.group(2)

        if metric == 'M':
            numG /= 1024
        elif metric == 'K':
            numG /= 1024**2
    else:
        memory_val = int(memory_val)
        numG = memory_val / (1024**3)

    return numG 


def runJob(cmd_job, collectl='collectl', interval=10, outfile=''):
    '''run job and use collectl to monitor system usage.
    collectl: location of collectl
    interval: interval in seconds to check system info.
    '''
    user_id = os.getuid()
    dat_file = outfile + '.collectl'
    print('temp file name', dat_file,'\n')
    cmd_job = cmd_job + '&>'+dat_file+'.log'
    cmd_collectl = f'{collectl} --procopts w  --procfilt u{user_id} --interval :{interval} --flush {interval} -sZ -oD > {dat_file}'
    worker_collectl = Worker(cmd_collectl)
    worker_collectl.start()
    time0 = time.time()
    worker_job = Worker(cmd_job)
    worker_job.start()

    worker_job.proc.wait()
    print('total running seconds', time.time() - time0)
    worker_collectl.stop()
    return dat_file



def parse_collectl_dat(collectl_dat_file, PROGS_OF_INTEREST):

    prev_time = 0
    prognames_dict = dict()

    """
    line formatting:
    0       #Date
    1       Time
    2       PID
    3       User
    4       PR
    5       PPID
    6       THRD
    7       S
    8       VSZ
    9       RSS
    10      CP
    11      SysT
    12      UsrT
    13      Pct
    14      AccuTime
    15      RKB
    16      WKB
    17      MajF
    18      MinF
    19      Command
    """
    ls_lines = []
    

    with open(collectl_dat_file) as f:
        for line in f:
            line = line.strip()
            vals = re.split("\s+", line,maxsplit=19)
            if len(vals) < 20:
                continue
            if not any([e in vals[-1] for e in PROGS_OF_INTEREST]):
                continue
            if line.startswith('#'):
                continue
            for progname in PROGS_OF_INTEREST:
                if progname in PROGS_OF_INTEREST:
                    break
            vals[1] = "{} {}".format(vals[0], vals[1])
            vals[0] = progname
            ls_lines.append(vals)

    df = pd.DataFrame(ls_lines)
    df.columns = 'progname Time PID  User     PR  PPID THRD S   VSZ   RSS CP  SysT  UsrT Pct  AccuTime  RKB  WKB MajF MinF Command'.split()
    df['Time'] = df['Time'].apply(lambda x:datetime.datetime.strptime(x, "%Y%m%d %H:%M:%S"))
    df['Pct'] = df['Pct'].astype(int)
    df['RSS'] = df['RSS'].apply(compute_GB)
    d = df.groupby(['progname','Time'])
    df_sum = pd.DataFrame()
    for (progname, t),tdf in d:
        df_sum.loc[t, progname+'__cpu'] = tdf['Pct'].sum()
        df_sum.loc[t, progname+'__memory'] = tdf['RSS'].sum()
    df_sum['running_time'] = (df_sum.index - df_sum.index.min()).seconds
    df_sum = df_sum.sort_values(by='running_time')
    df_sum = df_sum.reset_index(drop=True)
    return df_sum


def main(cmd_job, collectl='collectl', interval=10, programs_of_interest='python', outfile='job_running_summary'):
    '''run commands of cmd_job. output a dataframe to outfile, with summary for each of programs_of_interest
    '''
    collectl_dat_file = runJob(cmd_job, collectl=collectl, interval=interval,outfile=outfile)
    PROGS_OF_INTEREST = programs_of_interest.split(',')

    df_sum = parse_collectl_dat(collectl_dat_file, PROGS_OF_INTEREST)
    df_sum.to_csv(outfile+'.tsv', sep='\t',index=None)


 
####################
description="""
generate time matrix from collectl dat file. 
Modified from https://github.com/trinityrnaseq/trinityrnaseq/blob/master/trinity-plugins/COLLECTL/util/collectl_dat_to_time_matrix.py
"""
 
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=description,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("--cmd", type=str,  required=True, help="cmd command to run")
    parser.add_argument("--interval", type=int, default=10, help="seconds intervals to check system usage. default: 10 seconds")
    parser.add_argument("--collectl", type=str, default='collectl', help="location of collectl, default: collectl")

    parser.add_argument("--out_prefix", type=str, default="job_running_summary", help="prefix for output files. default job_running_summary")

    parser.add_argument("--programs", type=str, default="", help="programs or keywords in the commandline to monitor. default ''. If more than one programs to monitor, join by ','. e.g. 'Trinity,Butterfly,samtools,inchworm,QuantifyGraph,BubbleUpClustering,jellyfish,GraphFromFasta,bowtie2'")

    f = parser.parse_args()
    cmd_job = f.cmd
    interval = f.interval
    programs = f.programs
    collectl = f.collectl
    out_prefix = f.out_prefix

    main(cmd_job = cmd_job, collectl=collectl, interval=interval, programs_of_interest=programs, outfile=out_prefix)
