import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def group_list(alist):
    adict = {}
    for line in alist:
        sp = line.split()[0]
        #sp = cmn.lastName(fn).split('_')[0]
        try:
            adict[sp].append(line)
        except KeyError:
            adict[sp] = [line]
    return adict            

def merge_lines(alist, blist):
    #take the largest total
    gdict = {}
    for line in alist+blist:
        ref = line.split()[1]
        try:
            gdict[ref].append(line)
        except KeyError:
            gdict[ref] = [line]
    new = []
    for ref in gdict:
        lines = gdict[ref]
        if len(lines) == 1:
            new += lines
        else:
            line = max(lines, key=lambda x: int(x.strip().split()[-2]))
            new.append(line)
    return new            

wdir = sys.argv[1]

cmd = 'cat %s/mapped_reads_count/*' % wdir
lines = cmn.cmd2lines(cmd)

gdict = group_list(lines)    

ddir = '/project/biophysics/Nick_lab/mtang/archive/step2_bwa_mapping/mapped_reads_count'

for sp in gdict:
    lines = gdict[sp]
    dn = '%s/%s_cov.count' % (ddir, sp)
    if os.path.exists(dn):
        print('merging old data for %s' % dn)
        old_lines = cmn.file2lines(dn)
        lines = merge_lines(lines, old_lines)
    lines.append('')
    cmn.write_lines(lines, dn)        
    print('finish archiving %s' % sp)
