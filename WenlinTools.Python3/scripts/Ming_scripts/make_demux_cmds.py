import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def parse_index_mapping(fn):
    global wdir
    fi = '/work/biophysics/mtang/SNP_calling/scripts/data/adaptor/index_and_adaptor'
    adict = {}
    for line in cmn.file2lines(fi):
        index, barcode, seq = line.strip().split()
        adict[index] = barcode
    
    taken = set([])
    new = []
    for line in cmn.file2lines(fn):
        try:
            sp, index = line.strip().split()
        except:
            continue

        if not index.isdigit():
            continue
        
        #check for duplication
        if index in taken:
            print('Error! duplicated index detected for %s! please fix!' % wdir)
        else:
            taken.add(index)
        
        new.append('%s\t%s\n' % (sp, adict[index]))

    dn = 'lib_index.txt'
    cmn.write_file(''.join(new), dn)
    return dn
        

def tell_libs():
    libs = cmn.cmd2lines('ls *R*q')
    R1, R2, R3 = None, None, None
    for lib in libs:
        if 'R1' in lib:
            R1 = lib
        elif 'R2' in lib:
            R2 = lib
        elif 'R3' in lib:
            R3 = lib

    if R1 == None or R2 == None or R3 == None:
        print('can not tell orders for fastq:')
        print('\n'.join(libs))
        print('please fix!')
        sys.exit()
    return R1, R2, R3        
        

#---------------main--------------------
wdirs = cmn.cmd2lines('ls -d */| grep -v dm;')

fjobs = []
for wdir in wdirs:
    
    wdir = wdir.rstrip('/')
    #1. check if index_mapping.txt exists
    if not os.path.exists('%s/index_mapping.txt' % wdir):
        print('Error! no Index file found for %s! please fix!' % wdir)
        continue
    
    #2. process index mapping
    outdir = '%s_demux' % (wdir)
    cmn.mkdir(outdir)

    os.chdir(wdir)
    findex = parse_index_mapping('index_mapping.txt')
    R1, R2, R3 = tell_libs()
    
    cmds = ['cd %s' % outdir]
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/demux.py ../%s/%s ../%s/%s ../%s/%s 1 &' % (wdir, findex, wdir, R1, wdir, R2)
    cmds.append(cmd)
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/demux.py ../%s/%s ../%s/%s ../%s/%s 2 &' % (wdir, findex, wdir, R3, wdir, R2)
    cmds.append(cmd)
    cmds.append('\nwait\n')    

    os.chdir('..')

    #make job 
    fcmd = 'dm%s.cmd' % wdir
    cmn.write_lines(cmds, fcmd)
    fjob = 'dm%s.job' % wdir
    cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -p 256GB > %s' % (fcmd, fjob)
    cmn.run(cmd)

    fjobs.append(fjob)

print('please use the following commands to submit jobs for demux:')

for fn in fjobs:
    print('sbatch %s' % fn)
