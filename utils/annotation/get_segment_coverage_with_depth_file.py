# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 14:55:49 2018

@author: ATPs
"""
description= '''
input two file, file_partition and file_depth
file_partition is the one generate with Qian or Jing's code
file_depth is the depth file from samtools depth, output all position
by default, use 32 CPUs
write segment_coverage_100 file in current folder
'''
print(description)

import sys
from multiprocessing import Pool
import numpy as np

file_partiton = sys.argv[1]
file_depth = sys.argv[2]

ls_scf = []
fo = open(file_depth)
line = fo.readline()
scf, pos, depth = line.split()
ls_scf.append([scf,[]])
ls_scf[-1][1].append(int(depth))

for line in fo:
    scf, pos, depth = line.split()
    if scf == ls_scf[-1][0]:#the same scf
        ls_scf[-1][1].append(int(depth))
    else:#different scf
        ls_scf[-1][1] = np.array(ls_scf[-1][1],dtype=np.int16)#convert list of depth to numpy.array
        ls_scf.append([scf,[]])
        ls_scf[-1][1].append(int(depth))

ls_scf[-1][1] = np.array(ls_scf[-1][1],dtype=np.int16)#finish reading the file, convert the last one of depth
fo.close()


dc_scf=dict(ls_scf)

def get_depth(line):
    '''
    line is a line in file_partion, which looks like 
        scaffold169669_cov76:1-101      100
        scaffold169669_cov76:101-201    100
    return value that is ready to write, which is
        'part_id    part_len    reads    depth' 
        reads = int(depth/150)
    '''
    part_id, part_len = line.split()
    part_scf, part_se = part_id.split(':')
    part_start, part_end = part_se.split('-')
    part_start = int(part_start)-1
    part_end = int(part_end)-1
    part_len = int(part_len)
    np_depths = dc_scf[part_scf][range(part_start,part_end)]
    depth = int(np_depths.mean())
    reads = int(depth * part_len/100)
    return '{}\t{}\t{}\t{}\n'.format(part_id, part_len, reads, depth)


fo = open(file_partiton)
ls_lines = fo.readlines()
pool = Pool(32)
info_depths = pool.map(get_depth, ls_lines)
pool.close()

fout = open('segment_coverage_100','w')
for e in info_depths:
    fout.write(e)
fout.close()


