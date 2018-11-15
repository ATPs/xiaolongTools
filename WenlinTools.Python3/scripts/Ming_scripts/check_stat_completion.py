import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def group_list(alist):
    adict = {}
    for fn in alist:
        sp = cmn.lastName(fn).split('_')[0]
        try:
            adict[sp].append(fn)
        except KeyError:
            adict[sp] = [fn]
    return adict            


def check_NA(label):
    fn = label + '_stat.report'
    line = cmn.file2lines(fn)[-1]
    items = line.strip().split()
    #return isWarnning
    if len(items) != 10 or 'NA' in items:
        return True
    else:
        return False
fq_list = cmn.getid(sys.argv[1])
vcf_list = cmn.getid(sys.argv[2])
samdir_list = cmn.getid(sys.argv[3])
vcfCov_dir = cmn.getid(sys.argv[4])[0]

report_files = cmn.cmd2lines('ls *_stat.report')
finished_labels = set([cmn.lastName(each).replace('_stat.report', '')
                for each in report_files])

refresh = any([each=='-r' for each in sys.argv])

fq_groups = group_list(fq_list)# group by sp 
vcf_groups = group_list(vcf_list)

isGood = True

warnning_list = []

cmds = []
for sp in vcf_groups:
    vcf_fns = vcf_groups[sp]
    fq_fns = fq_groups[sp]

    #find the bwadir
    if len(samdir_list) == 1:
        samdir = samdir_list[0]
    else:
        #find the right dir
        isFound = False
        for samdir in samdir_list:
            checkFn = '%s/mapped_reads_count/%s_cov.count' % (samdir, sp)
            if os.path.exists(checkFn):
                isFound = True
                break
        if not isFound:
            print('Error! can not detect sam dir for %s' % sp)
            sys.exit()
    
    for vcf_fn in vcf_fns:
        label = cmn.lastName(vcf_fn).replace('_snp_step2.vcf', '')
        if label in finished_labels:
            isWarnning = check_NA(label)
            if isWarnning:
                warnning_list.append(label)
            continue
            
        isGood = False
        cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/generate_mapping_stat.py %s ' % ','.join(fq_fns)
        cmd += '%s %s %s ' % (vcf_fn, samdir, vcfCov_dir)
        cmds.append(cmd)

if isGood:
    if len(warnning_list) == 0:
        print('good news! everything looks good!')
    else:
        print('***************************************')
        print('Warnning: Please Email the following suspicious IDs to Wenlin:')
        print('\n'.join(warnning_list))
        print('***************************************')
        
else:
    cmds.append('')
    dn = 'statAdd.cmds'
    cmn.write_lines(cmds, dn)
    
    print('Error!!!!!')
    print('There are still %s report missing' % (len(cmds) - 1))
    print('please use following command to submit jobs')
    print('\n>>> /work/biophysics/mtang/SNP_calling/scripts/submit_jobs.py %s [#node] statAdd -p 256GB\n' % dn)
    print('-p specifies the partition it submitted to')
    print('[#node] is the number of nodes and should be adjusted according to number of lines in %s' % dn)
    print('\n[IMPORTANT]Please run this check again upon the job completion.')

