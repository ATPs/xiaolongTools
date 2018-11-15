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
import barcode_processing as bp


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    fns = cmn.cmd2lines('ls */sampleInfo.baits')

    splist = set([])
    for fn in fns:
        for line in cmn.file2lines(fn):
            a, b, c = line.strip().split()
            name = '%s_%s' % (a, b)
            splist.add(a)

    conn = bp.lock_database()
    for sp in splist:
        namelist = bp.select_barcodes("sampleID = '%s'" % sp, conn)
        for each in namelist:
            name = each[5]
            print(name)
            bp.delete_barcode_by_name(name, conn)

    bp.close_database(conn)
