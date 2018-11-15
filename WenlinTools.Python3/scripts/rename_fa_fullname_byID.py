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

    #t = ete3.Tree(cmn.txt_read(fn).replace('[&U]', ''))

    appear = {}
    table = []
    lines = cmn.file2lines(fn)
    info = []
    for line in lines:
        if line[0] =='>':
            sp = line[1:].split('_')[0]
            print(sp)
            try:
                newline = '>' + nameDict[sp].replace('\t', '_').replace(' ', '_')
            except:
                newline = line
            info.append(newline)
        else:
            info.append(line)

    info.append('')
    info = '\n'.join(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)
    cmn.write_file(''.join(table), dn + '.nameTable')



