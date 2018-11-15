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
import os
from fullname_lib import get_names_4barcode
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':

    #infoLines = cmn.cmd2lines('head -n 1 sampleRun*/rescued_read_assembled_mis1*.txt')
    IDlist = cmn.cmd2lines("ls -d sampleRun_* |grep -v fake|cut -d '_' -f 2")

    nameDict = get_names_4barcode()

    for ID in IDlist:
        items = nameDict[ID].replace('?', '').split()
        ID, genus, sp = items[:3]
        print('sampleInfo', ID, genus, sp)

        fn = 'sampleRun_%s/good_read_assembled.txt' % ID
        #label = '%s_%s' % (genus, sp)
        cmd = 'head %s -n 2| grep %s' % (fn, genus)
        print(cmd)
        info = cmn.cmd2info(cmd).strip()
        if info == '':
            print('please re-run', ID, genus, sp)


