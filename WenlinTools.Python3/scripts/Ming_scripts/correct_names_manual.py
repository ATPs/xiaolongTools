import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

vcf_list = [os.path.abspath(fn) for fn in cmn.getid(sys.argv[1])]

ref = 'Junonia_v2'

for fn in vcf_list:
    items = fn.split('/')
    sp = cmn.lastName(fn).split('_')[0]
    
    parent_dir = '/'.join(items[:-1])
    new_vcf = '%s/%s_%s_snp_step2.vcf' % (parent_dir, sp, ref)
    cmd = 'mv %s %s' % (fn, new_vcf)
    print(cmd)
    parent_dir = '/'.join(items[:-2])
    newdir = '%s/%s_%s' % (parent_dir, sp, ref)
    parent_dir = '/'.join(items[:-1])
    cmd = 'mv %s %s' % (parent_dir, newdir)
    print(cmd)
