import pandas as pd
import re
import time
import threading
import random
import os
import glob
import psutil
import shutil

class Worker(threading.Thread):
    def __init__(self, shell):
        self.shell = shell
        threading.Thread.__init__(self)
    def run(self):
        os.system(self.shell)
        time.sleep(random.randint(10, 100) / 1000.0)

def availableMemoryEnough(memory_min):
    '''
    check if available memory is larger than memory_min
    return True or False
    memory_min is float number in GB
    '''
    memory = psutil.virtual_memory()
    available = memory.available / (1024 ** 3)
    return available > memory_min


def multiThread(jobs, thread_max=1, memory_min=None, time_sleep=2):
    '''
    run jobs in multiple thread, with maximum thread_max jobs running, and only start a new job if the available memory is larger than memory_min (in GB)
    jobs is a folder or a file, each line is a job to run
    time_sleep is the gap to wait to check number of running jobs and available memory
    '''
    thread_count = thread_max
    if isinstance(jobs, list):
        jobs = jobs
    elif os.path.isdir(jobs):
        jobs = glob.glob(os.path.join(jobs,'*'))
    else:
        jobs = open(jobs).read().strip().split('\n')
    jobs =[e for e in jobs if e[0]!='#'] #if a line begin with #, do not run that line
    
    threading_dict = {}
    threading_count = 1
    for lines in jobs:
        threading_dict[threading_count] = Worker(lines)
        threading_count += 1

    for keys in threading_dict:
        threading_dict[keys].start()
        time.sleep(time_sleep)
        while True:
            if threading.active_count() > thread_count or not availableMemoryEnough(memory_min):
                time.sleep(2)
            else:
                break
    # wait all threads finish
    for keys in threading_dict:
        threading_dict[keys].join()
    print('all threads finished')


def getIntervals(file_interval):
    df_intervals = pd.read_csv(file_interval, sep='\t',header=None,comment='@')
    # keep only first 3 columns
    df_intervals = df_intervals[df_intervals.columns[:3]]
    df_intervals.columns=['chr','start','end']
    # sort by interval length, so that longer intervals will be processed first to save time
    df_intervals['len'] = df_intervals['end'] - df_intervals['start']
    df_intervals = df_intervals.sort_values(by=['len'], ascending=False)
    
    df_intervals['interval'] = df_intervals.apply(lambda x:'{}:{}-{}'.format(x['chr'], x['start'], x['end']), axis=1)
    
    return list(df_intervals['interval'])


def runHaplotypeCaller(file_interval, cmd, threads, memory_min, time_sleep=0):
    file_output = (re.findall('-O *.*? ',cmd) + re.findall('--output *.*? ',cmd))[0].split()[1]
    folder_output = os.path.dirname(file_output)
    folder_temp = file_output+'__HaplotypeCaller__temp'
    if not os.path.exists(folder_temp):
        os.makedirs(folder_temp)
    ls_intervals = getIntervals(file_interval)
    ls_cmds = []
    for interval in ls_intervals:
        t_out = os.path.join(folder_temp, interval+'.g.vcf.gz')
        t_log = os.path.join(folder_temp, interval+'.log')
        if not os.path.exists(t_out):
            ls_cmds.append(cmd.replace(file_output, t_out) + ' -L ' + interval +' &>' + t_log)
    #ls_cmds = [cmd.replace(file_output, t_out) + ' -L ' + interval +' &>' + t_log for interval in ls_intervals]
    
    multiThread(jobs=ls_cmds, thread_max=threads, memory_min=memory_min, time_sleep=time_sleep)
    
    # get all results g.vcf.gz
    files_gvcf = glob.glob(folder_temp+'/*.g.vcf.gz')
    file_gvcfs = os.path.join(folder_temp,'all_gvcfs.txt')
    open(file_gvcfs,'w').write('\n'.join(files_gvcf))
    # combine gatk results
    cmd_gatk = cmd.split('HaplotypeCaller')[0]
    cmd_gatk = cmd_gatk + f' GatherVcfs -I {file_gvcfs} -O {file_output}  -RI true --CREATE_INDEX true  &>{folder_temp}/GatherVcfs.log'
    
    run_status = os.system(cmd_gatk)
    if run_status == 0:
        shutil.rmtree(folder_temp)
        print('finished')
    else:
        print('something wrong, check tempfolder', folder_temp)




description = """
    run HaplotypeCaller in multiple thread, with maximum thread_max jobs running, and only start a new job if the available memory is larger than memory_min (in GB). below are example of input
file_interval = '/gpfs/gpfs/staging/jx76-003/xc/share/HumanGenome/GATK_resources/resources_broad_hg38_v0_wgs_calling_regions.hg38.interval_list'
cmd = '''gatk --java-options "-Dsamjdk.compression_level=5 -Xms5G -Djava.io.tmpdir=/home/scratch"  HaplotypeCaller      -R /gpfs/gpfs/staging/jx76-003/xc/share/HumanGenome/GRCh38DH_1KGP/GRCh38_full_analysis_set_plus_decoy_hla.fa      -I /gpfs/gpfs/staging/jx76-003/xc/20200812perGeno/practice/SRR1573206.GatherBamFiles.bam     -O /gpfs/gpfs/staging/jx76-003/xc/20200812perGeno/practice/SRR1573206.g.vcf.gz      -contamination 0      -G StandardAnnotation -G StandardHCAnnotation  -G AS_StandardAnnotation  -GQB 10 -GQB 20 -GQB 30 -GQB 40 -GQB 50 -GQB 60 -GQB 70 -GQB 80 -GQB 90    -ERC GVCF'''
threads=24
memory_min=30
time_sleep=0
"""

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--interval', help = 'a file of intervals to split the genome for -L options of HaplotypeCaller', required=True)
    parser.add_argument('-c','--cmd', help = 'the cmd scripts to run HaplotypeCaller', required=True)
    parser.add_argument('-t','--threads', help = 'max threads to use, default 12', type=int, default=12)
    parser.add_argument('-m','--min_mem', help = 'minimum memory in GB to start a new thread, default 15', default=15, type=float)
    parser.add_argument('-s','--time_sleep', help = 'time_sleep is the gap to wait to check number of running jobs and available memory, default 30 seconds', default=30, type=float)
    f = parser.parse_args()
    runHaplotypeCaller(file_interval=f.interval, cmd=f.cmd, threads=f.threads, memory_min=f.min_mem, time_sleep=f.time_sleep)
