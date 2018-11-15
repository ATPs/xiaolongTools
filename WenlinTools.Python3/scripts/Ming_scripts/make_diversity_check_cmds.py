import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def group_fastqs(alist):
    adict = {}
    for each in alist:
        key = cmn.lastName(each).split('_')[0]
        try:
            adict[key].append(each)
        except KeyError:
            adict[key] = [each]
    
    #check and sort
    keys = list(adict.keys())
    for key in keys:
        each = adict[key]
        if len(each) != 2:
            print('Error! number of libs is wrong for %s' % key)
            print('below are the detected libs:')
            print('\n'.join(each))
            print('Please fix!')
            sys.exit()
        each.sort()
        adict[key] = each
    return adict        

def group_old_fastqs(alist, current_batch):
    adict = {}
    for each in alist:
        key = cmn.lastName(each).split('_')[0]
        
        batchID = each.split('/')[-2]
        if batchID == current_batch:
            continue

        try:
            adict[key].append(each)
        except KeyError:
            adict[key] = [each]
    
    #check and sort
    keys = list(adict.keys())
    for key in keys:
        each = adict[key]
        newlist = [[], []]
        each.sort()
        for line in each:
            if 'R1' in line:
                newlist[0].append(line)
            elif 'R2' in line:
                newlist[1].append(line)
        adict[key] = newlist
    return adict        

def transfer_archive_libs(key, oR1s, oR2s):    
    outdir = os.path.abspath('tmp_transfer')
    cmn.mkdir(outdir)
    allFns = oR1s + oR2s

    dns = []
    for fn in allFns:
        dn = '%s/%s' % (outdir, cmn.lastName(fn))
        if not os.path.exists(dn) and (not os.path.exists(dn[:-3])):
            cmd = 'rsync butterfly@toxea.swmed.edu:%s %s/' % (fn, outdir)
            print('transfering file %s' % fn)
            cmn.run(cmd)
        dns.append(dn)    

    for fn in dns:
        if os.path.exists(fn):
            cmd = 'gunzip %s' % fn
            cmn.run(cmd)
    
    R1s = cmn.cmd2lines('ls %s/%s_R1*.fastq' % (outdir, key))
    R2s = cmn.cmd2lines('ls %s/%s_R2*.fastq' % (outdir, key))
    return R1s, R2s

def merge_libs(R1, R2, oR1s, oR2s):
    outdir = os.path.abspath('combined_raw_libs')
    cmn.mkdir(outdir)

    info = [('R1', R1, oR1s), ('R2', R2, oR2s)]
    combined = {}
    for label, newR, oRs in info:
        print('merging libs for %s' % cmn.lastName(newR))
        cmd = 'cp %s %s' % (newR, outdir)
        cmn.run(cmd)

        dn = '%s/%s' % (outdir, cmn.lastName(newR))
        combined[label] = dn

        for oldR in oRs:
            cmd = 'cat %s >> %s' % (oldR, dn)
            cmn.run(cmd)
    
    return combined['R1'], combined['R2']

#---------------main--------------------
fastqs = cmn.getid(sys.argv[1])
current_batch = fastqs[0].split('/')[-2]

fastq_dict = group_fastqs(fastqs)

#search for old fastq in butterfly account
old_libs = cmn.cmd2lines('ssh butterfly@toxea.swmed.edu "ls /archive/butterfly/raw_data/*/*fastq.gz"')
old_fastq_dict = group_old_fastqs(old_libs, current_batch)

cmds = []
for key in fastq_dict:
    R1, R2 = fastq_dict[key]
    try:
        oR1s, oR2s = old_fastq_dict[key]
        print('detected old libs, processing data for %s' % key) 
        
        oR1s, oR2s = transfer_archive_libs(key, oR1s, oR2s)
        print(oR1s, oR2s)

        R1, R2 = merge_libs(R1, R2, oR1s, oR2s) 
    except KeyError:
        pass
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/get_short_wl.py %s %s %s' % (R1, R2, key)
    cmds.append(cmd)

    #make job 
    fcmd = 'diversity.cmds'
    cmn.write_lines(cmds, fcmd)

cmn.mkdir('divr_reports')    
line = 'total_raw_data\tuniq_reads\ttotal_reads\tpercent\n'
cmn.write_file(line, 'divr_reports/0000000_header')
