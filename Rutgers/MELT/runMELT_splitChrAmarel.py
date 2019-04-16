# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 18:00:56 2019

@author: ATPs
"""

#file_bam = '/projects/genetics/GTEx/v7/sorted_bams/coverage_processed_bams/SRR3485450_hg38_sorted.bam'
temp_folder ='/scratch/xc278/MELT/'

import os
import sys

file_bam = sys.argv[1]

srr = os.path.basename(file_bam).split('_')[0]

#create folder
folder_work = temp_folder + srr
if os.path.exists(folder_work):
    os.system('rm -rf '+folder_work)

os.makedirs(folder_work)
os.system(f'samtools view -H {file_bam} >{folder_work}/chr.txt')

#get chromosomes longer than 100000
lines = open(f'{folder_work}/chr.txt').readlines()
lines = [e for e in lines if e.startswith('@SQ')]
dc_chr = {}
for e in lines:
    es = e.split()
    dc_chr[es[1][3:]] = int(es[2][3:])
dc = {k:v for k,v in dc_chr.items() if v> 1000000}

#split bam file
file_cmds = f'{folder_work}/cmds.txt'
ls_cmds = [f'samtools view -bh -o {folder_work}/{chromosome}.bam {file_bam} {chromosome}' for chromosome in dc]
open(file_cmds,'w').write('\n'.join(ls_cmds))
os.system(f'python3 /cache/home/xc278/w/GitHub/xiaolongTools/multiThread.py 24 {file_cmds}')

#combine small bam files
small_chrs = [f'{folder_work}/{chromosome}.bam' for chromosome in dc if len(chromosome) > 5 or chromosome == 'chrY']
os.system(f'samtools merge {folder_work}/other_ori.bam ' + ' '.join(small_chrs))
chromosomes = [k for k in dc if len(k) <= 5 and k != 'chrY']
chromosomes.append('other')
for f in small_chrs:
    os.remove(f)
os.system(f'samtools sort -o {folder_work}/other.bam -@ 24 {folder_work}/other_ori.bam')
os.remove(f'{folder_work}/other_ori.bam')

#index bam files
ls_cmds = [f'samtools index {folder_work}/{chromosome}.bam' for chromosome in chromosomes]
open(file_cmds,'w').write('\n'.join(ls_cmds))
os.system(f'python3 /cache/home/xc278/w/GitHub/xiaolongTools/multiThread.py 24 {file_cmds}')

# run MELT for each chromosome
ls_cmds = [f'/cache/home/xc278/p/java/jre1.8.0_201/bin/java -jar /cache/home/xc278/p/MELT/MELTv2.1.5/MELT.jar Genotype -h /scratch/xc278/MELT/GRCh38_full_analysis_set_plus_decoy_hla.fa -bamfile {folder_work}/{chromosome}.bam -t /cache/home/xc278/p/MELT/MELTv2.1.5/me_refs/Hg38/ALU_MELT.zip -p /scratch/xc278/MELT/MELT_ALU_combined -w {folder_work}/{chromosome}' for chromosome in chromosomes]
open(file_cmds,'w').write('\n'.join(ls_cmds))
os.system(f'python3 /cache/home/xc278/w/GitHub/xiaolongTools/multiThread.py 12 {file_cmds} 60')

#combine MELT result
files_result = [f'{folder_work}/{chromosome}/{chromosome}.ALU.tsv' for chromosome in chromosomes]
outfile = f'{folder_work}.ALU.tsv'
ls_results = [open(f) for f in files_result]
ls_good = []
for results in zip(*ls_results):
    good_results = [e for e in results if e.split()[1] != './.:-0,-0,-0']
    if len(good_results) == 0:
        good_result = results[0]
    else:
        good_result = good_results[0]
    ls_good.append(good_result)

open(outfile,'w').write(''.join(ls_good))
os.system('rm -rf '+ f'{folder_work}')



