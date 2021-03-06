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

fq_list = cmn.getid(sys.argv[1])
vcf_list = cmn.getid(sys.argv[2])
samdir_list = cmn.getid(sys.argv[3])
vcfCov_dir = cmn.getid(sys.argv[4])[0]

refresh = any([each=='-r' for each in sys.argv])

fq_groups = group_list(fq_list)# group by sp 
vcf_groups = group_list(vcf_list)

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
        cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/generate_mapping_stat.py %s ' % ','.join([eachone.strip() for eachone in fq_fns])
        cmd += '%s %s %s ' % (vcf_fn, samdir, vcfCov_dir)
        cmds.append(cmd)

cmds.append('')
dn = 'stat.cmds'
cmn.write_lines(cmds, dn)
