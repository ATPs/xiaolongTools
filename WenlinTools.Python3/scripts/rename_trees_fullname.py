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

from fullname_lib import get_names
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



def rename_tree(line):
    global nameDict
    t = ete3.Tree(line.replace('[&U]', ''))

    appear = {}
    table = []
    for node in t:
        name = node.name
        #sp = name.split('_')[0].split('.')[0]
        sp = name.split('_')[0]
        if sp not in appear:
            appear[sp] = 1
            table.append('%s\n' % ( nameDict[sp].replace(' ', '_')))
        else:
            appear[sp] += 1

        new_name = '%s_cp%s' % (nameDict[sp], appear[sp])
        #new_name = name.replace(sp, nameDict[sp])
        node.name = new_name

    info = t.write(format=1)
    return info


if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    nameDict = get_names()

    new = []
    for line in cmn.file2lines(fn):
        newTree = rename_tree(line)
        new.append(newTree)

    new.append('')
    dn = fn + '.renamed'
    cmn.write_lines(new, dn)



