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
import ete3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        ftree, faln=sys.argv[1:]
    except:
        print("Usage: *.py ftree fa", file=sys.stderr)
        sys.exit()



    tree = ete3.Tree(ftree)

    takenIDs = []
    with open(faln) as fp:
        for line in fp:
            if line[0] == '>':
                ID = line[1:].strip()
                takenIDs.append(ID)

    tree.prune(takenIDs)

    dn = cmn.lastName(ftree).replace('.tre', '') + '_prune.tre'
    cmn.write_file(tree.write(), dn)
