# -*- coding: utf-8 -*-
"""
Created on Sun Dec 23 18:24:38 2018

@author: ATPs
"""

file_genome = '/home/xcao/w/20181205Hermeuptychia/20181208QianPipeline/8_improve_assembly/wenlin3303min200filter.fa'
workfolder = '/home/xcao/w/20181205Hermeuptychia/20181208QianPipeline/9_identify_repeats_by_coverage'
file_bam = '/home/xcao/w/20181205Hermeuptychia/20181208QianPipeline/7_align_reads/merged.bam'

import os
from Bio import SeqIO
from multiprocessing import Pool

def get_scf_len(file_genome, outfile=None):
    '''
    given a file_genome, output a file with scaffold length
    for each line, with scf_id, scf_len, scf_len_noN
    outfile will file_genome+'.scf_size' if not provided
    '''
    if outfile is None:
        outfile = file_genome + '.scf_size'
    
    fout = open(outfile,'w')
    for s in SeqIO.parse(file_genome,'fasta'):
        seq = str(s.seq)
        l = len(seq)
        seq.upper()
        l_noN = l -seq.count('N')
        fout.write('%s\t%d\t%d\n'%(s.id, l, l_noN))
    
    fout.close()

def partition_scafs_step1(file_scfsize, file_genome, outfile):
    '''
    '''
    scfs = []
    fp = open(file_scfsize,"r").readlines()
    for line in fp:
        words = line.split()
        scaf = words[0]
        scfs.append(scaf)
    
    scafseqs = {}
    fp = open(file_genome,"r").read().split(">")
    for item in fp:
        if item:
            lines = item.split("\n")
            scaf = lines[0]
            seq = ""
            for line in lines[1:]:
                seq += line
            scafseqs[scaf] = seq
    
    fout = open(outfile,'w')
    
    for scf in scfs:
        if scf not in scafseqs:
            continue
        subseqs = scafseqs[scf].split("N")
        start = 0
        for subseq in subseqs:
            start += 1
            end = start + len(subseq)
            if end == start:
                continue
            elif end - start <= 100:
                fout.write(scf + ":" + str(start) + "-" + str(end) + "\t" + str(end - start) + '\n')
            else:
                substart = 0
                for i in range(len(subseq)//100):
                    substart = i*100
                    fout.write(scf + ":" + str(start + substart) + "-" + str(start + substart + 100) + "\t100" + '\n')
                fout.write(scf + ":" + str(start + substart + 100) + "-" + str(end) + "\t" + str(end - start - substart - 100) + '\n')
            start += len(subseq)
    
    fout.close()


def partition_scafs_step2(file_in, outfile):
    '''
    '''
    fout = open(outfile,'w')
    segments = []
    fp = open(file_in,"r").readlines()
    for line in fp:
        words = line.split()
        if int(words[1]) >= 100:
            segments.append(words)
        else:
            seginfo = words[0].split(":")
            segstart = seginfo[1].split("-")[0]
            segend = seginfo[1].split("-")[1]
            if segments[-1][0].split(":")[0] == seginfo[0]:
                segments[-1][0] = segments[-1][0].split("-")[0] + "-" + str(segend)
                segments[-1][1] = str(int(segments[-1][1]) + int(words[1]))
            else:
                segments.append(words)
    for segment in segments:
        fout.write(segment[0] + "\t" + segment[1] + '\n')
    fout.close()


def get_potenital_repeats(file_genome, file_bam, workfolder = '.', samtools='samtools', threads=32):
    '''
    given a genome sequence, get potential repeats in genome 
    '''
    file_scfsize = os.path.join(workfolder,'scf_size')
    file_partition_scfs_1 = os.path.join(workfolder, 'partition_scafs_1')
    file_partition_scfs_2 = os.path.join(workfolder, 'partition_scafs_2')
    get_scf_len(file_genome=file_genome, outfile=file_scfsize)
    partition_scafs_step1(file_scfsize=file_scfsize, file_genome=file_genome, outfile=file_partition_scfs_1)
    partition_scafs_step2(file_in=file_partition_scfs_1, outfile=file_partition_scfs_2)
    
    os.system('cd')
    commandline = '{samtools} depth -r {region}'
    
    
    
    
    