import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def filter_best_fastq(fns):
    #group them by ID
    gdict = {}
    for fn in fns:
        sp = cmn.lastName(fn).split('_')[0]
        try:
            gdict[sp].append(fn)
        except KeyError:
            gdict[sp] = [fn]
    #check how many different parent dict for each one
    newlist = []
    for sp in gdict:
        fns = gdict[sp]
        pdirs = {}
        for fn in fns:
            pdir = '/'.join(fn.split('/')[:-1])
            try:
                pdirs[pdir].append(fn)
            except KeyError:
                pdirs[pdir] = [fn]
        if len(pdirs) == 1:
            newlist += fns
        else:
            #if multiple data, then
            #1. check to take the one with the biggest file size
            maxFns = (0, None)
            for pdir in pdirs:
                subFns = pdirs[pdir]
                size = sum([cmn.filesize(each) for each in subFns])
                if size > maxFns[0]:
                    maxFns = (size, subFns)
            newlist += maxFns[1]
    return newlist            
            

wdir = sys.argv[1]

os.chdir(wdir)
#find the two input files

fns = cmn.cmd2lines('ls ')
auto_files = set(['mapping_info.txt', 'require_SNPs.dict.pk'])
fns = set(fns) - auto_files

fq = ''
fref = ''

for fn in fns:
    lines = cmn.file2lines(fn)
    line = lines[0]
    items = line.split()
    #print items
    if len(items) == 1:
        if 'fq' in line or ('fastq' in line):
            fq = fn
    elif len(items) >= 2:
        if 'fq' not in line and ('fastq' not in line) and ('fa' in line):
            fref = fn

if fq == '':
    print('Error! can not find fastq list file!')
    sys.exit()
else:
    print('guessing fastq file to be %s' % fq)

if fref == '':
    print('Error! can not find ref table file!')
    sys.exit()
else:
    print('guessing ref table file to be %s' % fref)


fq_all = '/project/biophysics/Nick_lab/mtang/archive/step1_info/fastq.filelist'
if os.path.exists(fq_all):
    aset = set(cmn.getid(fq_all))
else:
    aset = set([])

bset = set(cmn.getid(fq))
newset = aset | bset

newset = filter_best_fastq(newset)

cmn.write_lines(newset, fq_all)

fref_all = '/project/biophysics/Nick_lab/mtang/archive/step1_info/refTable.txt'
if os.path.exists(fref_all):
    aset = set(cmn.getid(fref_all))
else:
    aset = set([])

bset = set(cmn.getid(fref))
newset = aset | bset
cmn.write_lines(newset, fref_all)
   

