# -*- coding: utf-8 -*-
"""
Created on Sun Dec 16 18:30:52 2018

@author: ATPs
"""
import subprocess
import os
import time

def getNumberOfRunningJobs(user = 's185491', keyword=''):
    '''
    get the number of jobs running for user with the keyword
    if user == 'all', get result for all users
    '''
    if user == 'all':
        running_txt = subprocess.check_output('squeue',shell=True).decode('utf-8')
    else:
        running_txt = subprocess.check_output('squeue -u '+user,shell=True).decode('utf-8')
    
    running_list = running_txt.strip().split('\n')
    running_list = running_list[1:]
    running_list = [e for e in running_list if keyword in e]
    return len(running_list)

def submitJobs(jobs,N=1,user = 's185491', keyword='',timesleep=60):
    '''
    jobs is a list of sbatch filenames to run, or a file with the location of these jobs to submit
    N is the number of jobs to run at the same time
    timesleep is the gap to check if N jobs were running. default, 60 seconds
    keyword is the keyword to check when use 'squeue'
    user is the user to check
    '''
    if isinstance(jobs, str):
        jobs = open(jobs).readlines()
        jobs = [e.strip('\n') for e in jobs]
    
    jobs = [e for e in jobs if not e.startswith('#')]
    print(len(jobs),'jobs to submit. jobs running at the same time with keyword:[',keyword,']were set to', N)
    
    
    for job in jobs:
        while True:
            jobsRunning = getNumberOfRunningJobs(user=user,keyword=keyword)
            if jobsRunning >= N:
                time.sleep(timesleep)
            else:
                break
        os.system('qsub '+job)
        print('submit job',job)
    
    print('done')

description = '''
    jobs is a list of sbatch filenames to run, or a file with the location of these jobs to submit
    lines start with '#' were ignored
    N is the number of jobs to run at the same time. default 1
    timesleep is the gap to check if N jobs were running. default, 60 seconds
    keyword is the keyword to check when use 'squeue'
    user is the user to check
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'input jobs file', required=True)
    parser.add_argument('-u','--user',help = 'username to check job with squeue, default s185491', default = 's185491')
    parser.add_argument('-N','--number', help = 'number of jobs running at the same time, default 1', default = 1, type=int)
    parser.add_argument('-k','--keyword', help = 'keyword to filter the jobs running. default, no filter', default = '')
    parser.add_argument('-t','--timesleep', help = 'check job status with the gap of timesleep seconds, default 60', default = 60, type=int)
    f = parser.parse_args()
    submitJobs(jobs=f.input, N=f.number,user = f.user, keyword=f.keyword,timesleep=f.timesleep)
