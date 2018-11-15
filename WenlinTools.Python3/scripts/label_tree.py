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
    ftree = sys.argv[1]#'Junonia_2017-07-10_ExaML_scafBT_cc.tre'
    fpick = sys.argv[2]#'filtered_table_step2.txt'
    pickedIDs = set([line.split()[0].split('_')[0] for line in cmn.file2lines(fpick)])
    print(pickedIDs)

    info = cmn.txt_read(ftree)

    for ID in pickedIDs:
        info = info.replace('%s_' % ID, '***%s_' % ID)


    dn = cmn.lastName(ftree) + '_marked.tre'
    cmn.write_file(info, dn)

