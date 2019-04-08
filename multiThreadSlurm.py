# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 14:57:32 2018

@author: ATPs
"""
import re
import time
import sys
import threading
import time
import random
import os
import glob


class Worker(threading.Thread):
    def __init__(self, shell):
        self.shell = shell
        threading.Thread.__init__(self)
    def run(self):
        os.system(self.shell)
        time.sleep(random.randint(10, 100) / 1000.0)
            
def main(Args):
    if not Args:
        print ("Usage: python <thread_count> <in_file>")
        sys.exit(0)
    thread_count = int(Args[0])
    in_file = Args[1]

    threading_dict = {}
    threading_count = 1
    for lines in open(in_file, "r"):
        if re.match("#", lines) or not re.search("[0-9a-zA-Z]+", lines):
            continue
        threading_dict[threading_count] = Worker(lines.rstrip())
        threading_count += 1

    for keys in threading_dict:
        threading_dict[keys].start()
        time.sleep(1)
        while True:
            if threading.active_count() > thread_count:
                time.sleep(2)
            else:
                break

def getJobs(jobs):
    '''
    jobs is a file, each line is a job that can run directly, or a list of commands, or a folder with jobs to run
    return a list of cmds to run
    '''
    if isinstance(jobs, list):
        print('input is a list of commandlines')
        return jobs
    if os.path.isfile(jobs):
        print('input is a filename')
        jobs = open(jobs).read()
        jobs = jobs.strip()
        jobs = jobs.split('\n')
        return jobs
    if os.path.isdir(jobs):
        print('input is a folder with scripts to run')
        jobs = glob.glob(jobs + '/*')
        jobs = ['bash '+e for e in jobs]
        return jobs
    print('something wrong with the input jobs')
    exit(-1)

def getNodes(nodes):
    '''
    nodes is a file, each line is a node name that can be accessed via ssh without a password, or a str of nodes description separated by ',', like elc-109;elc-110, or node074;node084; default None, work in current node
    return None, or a list of node names
    '''
    if nodes is None:
        print('nothing provided, input nodes is None, run in current node')
        return None
    if os.path.isfile(nodes):
        print(nodes, ':input nodes is a file')
        nodes = open(nodes).read().splitf()
        return nodes
    else:
        print(nodes, ': input nodes is a str of nodes')
        nodes = nodes.split(',')
        return nodes

def getWorkNode(ls_jobs, nodes, threads):
    '''
    decide which node to run the work_job
    a helper function
    '''
    active_jobs = [e for e in ls_jobs if e.isAlive()]
    dcNodes = {n:0 for n in nodes}
    for job in active_jobs:
        job_node = job.getName()
        dcNodes[job_node] += 1
    for k,v in dcNodes.items():
        if v < threads:
            return k
    return None

def submitJobs(jobs, nodes=None, threads=1, time_sleep=60):
    '''
    run jobs in different nodes so that each nodes have threads jobs running
    '''
    print('jobs',jobs)
    print('nodes',nodes)
    jobs = getJobs(jobs)
    nodes = getNodes(nodes)
    print('jobs',jobs)
    print('nodes',nodes)
    
    ls_jobs = [Worker(job) for job in jobs]
    
    if nodes is None:
        for work_job in ls_jobs:
            work_job.start()
            time.sleep(time_sleep)
            while True:
                if threading.active_count() > threads:
                    time.sleep(2)
                else:
                    break
    else:
        for work_job in ls_jobs:
            while True:
                work_node = getWorkNode(ls_jobs, nodes, threads)
                print(work_node)
                if work_node is None:
                    print('sleep',time_sleep)
                    time.sleep(time_sleep)
                else:
                    work_job.setName(work_node)
                    work_job.shell = """ssh {work_node} '{work_shell}'""".format(work_node=work_node, work_shell=work_job.shell)
                    print(work_job.shell)
                    work_job.start()
                    break
                
            

description = '''
    input is a file, each line is a job that can run directly, or a list of commands, or a folder with jobs to run
    nodes is a file, each line is a node name that can be accessed via ssh without a password, or a str of nodes description separated by ',', like elc-109;elc-110, or node074;node084; default None, work in current node
    t is the number of jobs to run in each node
    s sleep time about how often to check the running status
    '''
if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-i','--input', help = 'iinput is a file, each line is a job that can run directly', required=True)
    parser.add_argument('-n','--nodes',help = '''a file, each line is a node name that can be accessed via ssh without a password, or a str of nodes description separated by ',', like elc-109;elc-110, or node074;node084; default None, work in current node''', default = None)
    parser.add_argument('-t','--threads', help = 'number of jobs to run in each node, default 1', default = 1, type=int)
    parser.add_argument('-s','--sleep', help = 'sleep time about how often to check the running status, default 60', default = 60, type=int)
    f = parser.parse_args()
    submitJobs(jobs=f.input, nodes=f.nodes, threads=f.threads,time_sleep=f.sleep)