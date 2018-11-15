import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

fromDir = os.path.abspath(sys.argv[1])
toDir = os.path.abspath(sys.argv[2])

wdirs = cmn.cmd2lines('ls %s | grep ^sampleRun_' % fromDir)

#toKeep = ['*.txt', '*.report', 'barcode_count', '*_contig.fa', 'denovo_barcode.fa', 'bait0_denovo.br']
toKeep = ['*.txt', '*.report', '*_contig.fa', 'denovo_barcode.fa', 'bait0_denovo.br']
for wdir in wdirs:
    eachToDir = '%s/%s' % (toDir, wdir)
    cmn.mkdir(eachToDir)

    eachFromDir = '%s/%s' % (fromDir, wdir)

    for fn in toKeep:
        cmd = 'cp %s/%s %s' % (eachFromDir, fn, eachToDir)
        print(cmd)
        cmn.run(cmd)





