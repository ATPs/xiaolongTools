# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 17:16:25 2019

@author: ATPs
"""

folder = '/fastio/xc278/2019Bee/maker2/hornet/step1/' #folder of oringal maker
threads = 24 # number of threads

pre_cmd = 'source activate /gpfs/gpfs/scratch/xc278/maker2 ' # commond to run before each individual commands
pyMultiple = '/home1/xc278/w/GitHub/xiaolongTools/multiThreadSlurm.py'

import os
from Bio import SeqIO

#get setting files
file_maker_exe = os.path.join(folder,'maker_exe.ctl')
file_maker_opts = os.path.join(folder, 'maker_opts.ctl')
file_maker_bopts = os.path.join(folder,'maker_bopts.ctl')

#generate workfolder
workfolder = os.path.join(folder,'maker_split_run')
if not os.path.exists(workfolder):
    os.makedirs(workfolder)

#get genome location
for line in open(file_maker_opts):
    if line.startswith('genome='):
        break
file_genome = line[7:].split('#')[0].strip()

#split genome sequences to parts so that each parts contains is longer than 100kb
#remove contigs shorter than 200bp
seqfolder = os.path.join(workfolder,'seq')
if not os.path.exists(seqfolder):
    os.makedirs(seqfolder)

n = 0
l = 0
outfolder = os.path.join(seqfolder,str(n))
os.makedirs(outfolder)
outfile = open(os.path.join(outfolder,str(n)),'w')
for s in SeqIO.parse(file_genome,'fasta'):
    slen = len(s.seq)
    if slen <=200:
        continue
    l += slen
    outfile.write('>'+s.id+'\n'+str(s.seq)+'\n')
    if l >= 100000:
        outfile.close()
        n += 1
        l = 0
        outfolder = os.path.join(seqfolder,str(n))
        os.makedirs(outfolder)
        outfile = open(os.path.join(outfolder,str(n)),'w')
outfile.close()

#prepare to run maker for each part
total_seg = n+1
commands = os.path.join(workfolder,'cmds')
fout = open(commands,'w')
for n in range(total_seg):
    for f in [file_maker_exe, file_maker_opts,file_maker_bopts]:
        outfolder = os.path.join(seqfolder,str(n))
        genome = os.path.join(outfolder,str(n))
        os.system(f'cp {f} {outfolder}/')
    fout.write(pre_cmd + '  &&  ' + f'cd {outfolder} && maker -genome {genome} -cpus 1 \n')
fout.close()

#run jobs
cmd = f'python {pyMultiple} -t {threads} -i {commands} -s 10'
print(cmd)
#os.system(cmd)