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
        fn=sys.argv[1]
        fname = sys.argv[2]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    nameDict = cmn.pickle_read(fname)

    t = ete3.Tree(fn)

    for node in t:
        name = node.name
        node.name = nameDict[name]


    dn = cmn.lastName(fn) + '.mapnamed'
    cmn.write_file(t.write(), dn)





