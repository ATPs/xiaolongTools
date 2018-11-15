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


    fIDs = '/project/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/finished_sampleIDs.txt'

    IDlist = set([line.split()[0].replace('[new]','') for line in cmn.file2lines(fIDs)])

    new = []
    for fq in cmn.file2lines(fn):
        sample = cmn.lastName(fq).split('_')[0]
        if sample in IDlist:
            print('removed finished: %s' % fq)
        else:
            new.append(fq)

    dn = fn + '.filtered'
    cmn.write_lines(new, dn)
