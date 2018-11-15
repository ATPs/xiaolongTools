import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

wdir = os.path.abspath(sys.argv[1].rstrip('/'))

fvcfs = cmn.cmd2lines('ls %s/*/*.vcf' % wdir)

refdir = '/work/biophysics/mtang/SNP_calling/indexed_references'

badones = []
for fvcf in fvcfs:
    label = fvcf.split('/')[-2]
    reflabel = '_'.join(label.split('_')[1:])
    finfo = '%s/%s_scafLength.txt' % (refdir, reflabel)
    if not os.path.exists(finfo):
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/assembly_scaf_length.py %s/%s.fa ' % (refdir, reflabel)
        cmn.run(cmd)
    
    infoline = cmn.cmd2info('tail -n 1 %s' % finfo).strip()
    Cscaf, Cindex = infoline.split()[:2]

    checkline = cmn.cmd2info('tail -n 1 %s' % fvcf).strip()
    scaf, index = checkline.split()[:2]

    if scaf != Cscaf or Cindex != index:
        print('Error! problematic vcf file for %s' % label)
        badones.append(label)


dn = 'bad_vcf.list'
cmn.write_lines(badones, dn)


