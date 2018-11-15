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
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    new = []
    count = 0
    for line in cmn.file2lines(fn):
        count += 1
        t = ete3.Tree(line.replace('[&U]', ''))
        for leaf in t:
            leaf.dist -= 0.45
            if leaf.dist < 0:
                addD = 0.1 - leaf.dist
                leaf.dist = 0.1

                parentNode = leaf.up
                try:
                    isChange = parentNode.isChange
                except:
                    isChange = False

                if not isChange:
                    parentNode.dist -= addD
                    parentNode.isChange = True



        info = t.write(format=1)
        new.append(info)

    new.append('')
    cmn.write_lines(new, fn + '.cut')



