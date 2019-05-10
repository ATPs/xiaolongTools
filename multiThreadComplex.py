# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:57:32 2018

@author: ATPs
"""

#! /usr/bin/python
# -*- coding: utf-8 -*-
###################################################
# Program: **.py
# Function:
# Version: V1.0
# Date: Fri Oct 25 11:30:50 2013
###################################################
import re
import time
import sys
import threading
import time
import random
import os
import glob
import psutil

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
    if os.path.isdir(jobs):
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

description = '''
    run jobs in multiple thread, with maximum thread_max jobs running, and only start a new job if the available memory is larger than memory_min (in GB)
    jobs is a folder or a file, each line is a job to run
    time_sleep is the gap to wait to check number of running jobs and available memory, default 30 seconds
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'a file of jobs to run or a folder with job files', required=True)
    parser.add_argument('-t','--threads', help = 'max threads to use, default 12', type=int, default=12)
    parser.add_argument('-m','--min_mem', help = 'minimum memory in GB to start a new thread, default 15', default=15, type=float)
    parser.add_argument('-s','--time_sleep', help = 'time_sleep is the gap to wait to check number of running jobs and available memory, default 30 seconds', default=30, type=float)
    f = parser.parse_args()
    multiThread(jobs=f.input, thread_max=f.threads, memory_min=f.min_mem,time_sleep=f.time_sleep)
