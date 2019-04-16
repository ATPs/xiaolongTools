# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 12:37:01 2019

@author: ATPs
"""

# based on the author of MELT, it is better to split the pre_geno file
#file_bam = '/projects/genetics/GTEx/v7/sorted_bams/coverage_processed_bams/SRR3485450_hg38_sorted.bam'
#temp_folder ='/scratch/xc278/MELT/'
#folder_combined = '/projectsp/f_jx76_1/GTEx_v6_v7_jointMELT/MELT_Alu_combined'
#genome = '/scratch/xc278/MELT/GRCh38_full_analysis_set_plus_decoy_hla.fa'
#thread=13

import os
import sys
import psutil
import glob
import numpy as np



def run_MELT_SplitGene(file_bam, temp_folder, folder_combined, genome, thread=12):
    '''
    run MELT GenoType for file_bam in temp_folder. 
    folder_combined is the output folder of the previous step
    genome is the reference genome file
    thread is the number of thread to use
    '''
    srr = os.path.basename(file_bam).split('_')[0]

    #create folder
    folder_work = temp_folder + srr
    if os.path.exists(folder_work):
        os.system('rm -rf '+folder_work)
    os.makedirs(folder_work)
    
    #cp bam file to temp_folder
    basename_bam = os.path.basename(file_bam)
    os.system(f'cp {file_bam} {folder_work}/{basename_bam}')
    os.system(f'cp {file_bam}.bai {folder_work}/{basename_bam}.bai')
    file_bam = f'{folder_work}/{basename_bam}'
    
    #get the *.pre_geno.tsv file
    file_PreGeno = glob.glob(folder_combined+'/*.pre_geno.tsv')[0]
    
    #split the *.pre_geno.tsv file to thread parts
    ls_ref = open(file_PreGeno).readlines()
    ls_num = np.array_split(range(len(ls_ref)), thread)
    for i in range(thread):
        folder_ref = folder_work+'/ref'+str(i)
        os.makedirs(folder_ref)
        file_ref = folder_ref + '/' + os.path.basename(file_PreGeno)
        open(file_ref,'w').write(''.join([ls_ref[e] for e in ls_num[i]]))
    
    #run MELT for each section
    JAVA = '~/p/java/jre1.8.0_201/bin/java'
    MELT = '~/p/MELT/MELTv2.1.5/MELT.jar'
    ALU_MELT = '~/p/MELT/MELTv2.1.5/me_refs/Hg38/ALU_MELT.zip'
    ls_cmds = [f'{JAVA} -jar {MELT} Genotype -h {genome} -bamfile {file_bam} -t {ALU_MELT} -p {folder_work}/ref{i} -w {folder_work}/work{i}' for i in range(thread)]
    file_cmds = f'{folder_work}/cmds.txt'
    open(file_cmds,'w').write('\n'.join(ls_cmds))
    os.system(f'python3 ~/w/GitHub/xiaolongTools/multiThread.py {thread} {file_cmds} 60')
        
    #combine the result
    folders_result = [f'{folder_work}/work{i}' for i in range(thread)]
    files_result = [glob.glob(f+'/*.tsv')[0] for f in folders_result]
    outfile = temp_folder + '/' + os.path.basename(files_result[0])
    open(outfile,'w').write(''.join(open(f).read() for f in files_result))
    os.system('rm -rf '+ f'{folder_work}')
    
description = '''
    run MELT GenoType for file_bam in temp_folder. 
    folder_combined is the output folder of the previous step
    genome is the reference genome file
    thread is the number of thread to use
'''

if __name__ == '__main__':
    import argparse
    print(description)
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-b','--file_bam', help = 'bam file location', required=True)
    parser.add_argument('-w','--temp_folder', help = 'temp folder to work', required=True)
    parser.add_argument('-c','--folder_combined', help = 'the output folder of the previous step', required=True)
    parser.add_argument('-g','--genome', help = 'genome location. genome need to be indexed', required=True)
    parser.add_argument('-t','--thread', help = 'inumber of CPUs to use, default 12', type=int, default=12)
    f = parser.parse_args()
    run_MELT_SplitGene(file_bam=f.file_bam, temp_folder=f.temp_folder, folder_combined=f.folder_combined, genome=f.genome, thread=f.thread)