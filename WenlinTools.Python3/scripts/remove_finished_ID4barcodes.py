#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    f_IDs = cmn.cmd2lines('ls /project/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/*_sampleIDs.txt')

    IDlist = set([])
    for each in f_IDs:
        IDlist = IDlist | set([line.split()[0].replace('[new]','') for line in cmn.file2lines(each)])


    processingIDs = []
    fcompare = cmn.cmd2lines('ls -d /project/biophysics/Nick_lab/wli/sequencing/BWA_barcodes/*/compare.check')
    fcompare += cmn.cmd2lines('ls -d /project/biophysics/Nick_lab/wli/sequencing/BWA_barcodes/*/*/compare.check')
    fcompare += cmn.cmd2lines('ls /project/biophysics/Nick_lab/mtang/BWA_barcodes/*/compare.check')

    for feach in fcompare:
        processingIDs += cmn.cmd2lines('cut -f 1 %s| grep -v sample' % feach)
    processingIDs = set(processingIDs)

    new = []
    for sample in cmn.file2lines(fn):
        #sample = cmn.lastName(fq).split('_')[0]
        if sample in IDlist:
            print('removed finished: %s' % sample)
        else:
            if sample in processingIDs:
                print('removed processing ID: %s' % sample)
            else:
                new.append(sample)

    dn = fn + '.filtered'
    cmn.write_lines(new, dn)
