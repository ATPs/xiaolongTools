import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def backup_vcf_coverage(wdir):
    ddir = '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/check_vcf_coverage'
    fns = cmn.cmd2lines('ls %s/*_vcf.cov' % wdir)
    
    #1. only back up the new version of cov file
    for fn in fns:
        print('processing %s...' % fn)
        lines = cmn.file2lines(fn)
        items = lines[-1].strip().split()
        if len(items) != 6:
            print('skip old format file %s' % fn)
            continue
        fnlabel = cmn.lastName(fn)
        dn = '%s/%s' % (ddir, fnlabel)
        if os.path.exists(dn):
            print('merging new and old data for %s' % fnlabel)
            covOld = float(cmn.file2lines(dn)[-1].split()[-2])
            cov = float(lines[-1].split()[-2])
            if cov > covOld:
                cmn.run('cp %s %s' % (fn, dn))
        else:
            cmn.run('cp %s %s' % (fn, dn))
    

def count_map_nonGap(fn):
    gap_chars = set(['N', '-', 'X'])
    count = 0
    with open(fn) as fp:
        for line in fp:
            if line.strip().split()[0] not in gap_chars:
                count += 1
    return count


def backup_map(wdir):
    ddir = '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/vcf2map'
    fns = cmn.cmd2lines('ls %s/*.map' % wdir)
    #_snp_step2_MITO.map
    for fn in fns:
        print('processing %s...' % fn)
        fnlabel = cmn.lastName(fn)
        #don't back up the ones without species
        items = fnlabel.replace('_snp_step2_MITO.map','').replace('_snp_step2.map', '').split('_')
        if len(items) == 1:
            print('skip the map without sp for %s' % fn)
            continue
        
        #get the least gapped one
        dn = '%s/%s' % (ddir, fnlabel)
        if os.path.exists(dn):
            print('merging new and old data for %s' % fnlabel)
            Nold = count_map_nonGap(dn)
            Nnew = count_map_nonGap(fn)
            if Nnew > Nold:
                cmn.run('cp %s %s' % (fn, dn))
        else:
            cmn.run('cp %s %s' % (fn, dn))
        

def count_fasta_nonGap(fn):
    gap_chars = set(['N', '-', 'X'])
    count = 0
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '>':
                continue
            
            N = len(line)
            for char in gap_chars:
                N = N - line.count(char)
            
            count += N
            #for char in line:
            #    if char not in gap_chars:
            #        count += 1
    return count                    



def backup_fasta(wdir):
    ddir = '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/map2fasta'
    fns = cmn.cmd2lines('ls %s/*_m2s.fa| grep -v all_genome' % wdir)
    for fn in fns:
        print('processing %s...' % fn)
        fnlabel = cmn.lastName(fn)
        #don't back up the ones without species
        items = fnlabel.replace('_snp_step2_MITO_m2s.fa','').replace('_snp_step2_m2s.fa', '').split('_')
        if len(items) == 1:
            print('skip the fasta without sp for %s' % fn)
            continue
        
        #get the least gapped one
        dn = '%s/%s' % (ddir, fnlabel)
        if os.path.exists(dn):
            print('merging new and old data for %s' % fnlabel)
            Nold = count_fasta_nonGap(dn)
            Nnew = count_fasta_nonGap(fn)
            if Nnew > Nold:
                cmn.run('cp %s %s' % (fn, dn))
        else:
            cmn.run('cp %s %s' % (fn, dn))


def count_final_stat(fn):
    line = cmn.file2lines(fn)[-1]
    N = line.count('NA')
    Ndata = float(line.strip().split()[4])#expected coverage
    return N, Ndata


def backup_finalStat(wdir):
    ddir = '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/final_stats/'
    fns = cmn.cmd2lines('ls %s/*.report| grep -v all_genome' % wdir)
    for fn in fns:
        print('processing %s...' % fn)
        fnlabel = cmn.lastName(fn)
        #don't back up the ones without species
        items = fnlabel.replace('_stat.report','').split('_')
        if len(items) == 1:
            print('skip the fasta without sp for %s' % fn)
            continue
        
        #get the one with least NA and more data amount
        dn = '%s/%s' % (ddir, fnlabel)
        if os.path.exists(dn):
            print('merging new and old data for %s' % fnlabel)
            Nold_na, Nold_data = count_final_stat(dn)
            Nnew_na, Nnew_data = count_final_stat(fn)
            if Nnew_na < Nold_na: #less NA
                cmn.run('cp %s %s' % (fn, dn))
            else:
                if Nnew_na == Nold_na:#same NA number
                    if Nnew_data > Nold_data:
                        cmn.run('cp %s %s' % (fn, dn))
                        
        else:
            cmn.run('cp %s %s' % (fn, dn))
    cmn.run('cd %s; cat *.report > allstat.txt' % ddir)            
        
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wdir = os.path.abspath(sys.argv[1])

print('processing coverge files...')
subdir = '%s/check_vcf_coverage' % wdir
if os.path.exists(subdir):
    backup_vcf_coverage(subdir)
else:
    print('Error! can not find %s' % subdir)


print('processing fasta files...')
subdir = '%s/map2fasta' % wdir
if os.path.exists(subdir):
    backup_fasta(subdir)
else:
    print('Error! can not find %s' % subdir)

print('processing map files...')
subdir = '%s/vcf2map' % wdir
if os.path.exists(subdir):
    backup_map(subdir)
else:
    print('Error! can not find %s' % subdir)

print('processing final report files...')
subdir = '%s/final_stats' % wdir
if os.path.exists(subdir):
    backup_finalStat(subdir)
else:
    print('Error! can not find %s' % subdir)

