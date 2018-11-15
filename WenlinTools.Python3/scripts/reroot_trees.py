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

    nodes = ['%s_cp1' % each
            for each in ['6278','6277','6276']]

    new = []
    for line in cmn.file2lines(fn):
        t = ete3.Tree(line.replace('[&U]', ''))
        outgroup = t.get_common_ancestor(nodes)
        t.set_outgroup(outgroup)

        info = t.write(format=1)
        new.append(info)

    new.append('')
    cmn.write_lines(new, fn + '.reroot')



