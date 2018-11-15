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



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    nameDict = get_names()
    #print nameDict['DZ30.949']

    t = ete3.Tree(cmn.txt_read(fn).replace('[&U]', ''))

    appear = {}
    table = []
    for node in t:
        name = node.name
        sp = name.split('_')[0].replace('flt', '').replace('hc', '').strip('\'"').replace('.fa', '').split('.Le')[0]
        print('spCheck', sp)
        try:
            fullname = nameDict[sp]
        except KeyError:
            fullname = sp
        if sp not in appear:
            appear[sp] = 1
            table.append('%s\t%s\n' % (sp, fullname))
        else:
            appear[sp] += 1

        new_name = '%s_cp%s' % (fullname, appear[sp])
        #new_name = name.replace(sp, nameDict[sp])
        node.name = new_name

    info = t.write()
    print(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)
    cmn.write_file(''.join(table), dn + '.nameTable')



