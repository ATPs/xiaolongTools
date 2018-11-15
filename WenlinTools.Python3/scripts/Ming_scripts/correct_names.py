import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

vcf_list = [os.path.abspath(fn) for fn in cmn.getid(sys.argv[1])]

for fn in vcf_list:
    items = fn.split('/')
    parent_dir = '/'.join(items[:-3])
    step2_dir = '%s/step2_bwa_mapping/mapped_reads_count' % parent_dir
    sp = cmn.lastName(fn).split('_')[0]
    lines = cmn.cmd2lines('grep %s %s/*' % (sp, step2_dir))
    maxRef = (None, 0)
    if len(lines) == 0:
        print('Error for %s' % fn)

    for line in lines:
        a, ref, mapN, totalN = line.strip().split()
        if int(mapN) > maxRef[1]:
            maxRef = [ref, int(mapN)]
    
    ref = maxRef[0]
    #ref = 'Junonia_v2_withMito'
    parent_dir = '/'.join(items[:-1])
    new_vcf = '%s/%s_%s_snp_step2.vcf' % (parent_dir, sp, ref)
    cmd = 'mv %s %s' % (fn, new_vcf)
    print(cmd)
    parent_dir = '/'.join(items[:-2])
    newdir = '%s/%s_%s' % (parent_dir, sp, ref)
    parent_dir = '/'.join(items[:-1])
    cmd = 'mv %s %s' % (parent_dir, newdir)
    print(cmd)
