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
            name = leaf.name
            dist = leaf.dist
            pName = '_'.join(name.split('_')[:-1])
            pNode = leaf.up
            try:
                isChange = pNode.isChange
            except:
                isChange = False

            if not isChange:
                pNode.name = pName
                pNode.dist += dist
                pNode.isChange = True
            leaf.detach()


        info = t.write(format=2)
        new.append(info)

    new.append('')
    cmn.write_lines(new, fn + '.colps')



