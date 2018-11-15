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
    #print nameDict['DZ30.949']

    rdict = {}
    order_list = []
    for line in cmn.file2lines(fn):
        name, seq = line.strip().split()
        leftI = line.find(seq)
        sp = name.split('_')[0].split(']')[-1]
        try:
            fullname = nameDict[sp].replace(' ', '_').replace('\t', '_')
        except KeyError:
            fullname = name
            print('WARNNING: find no name for %s' % sp)
        order_list.append(fullname)
        rdict[fullname] = seq

    info = []
    nameL = max(list(map(len,list(rdict.keys())))) + 4
    for name in order_list:
        line = ('{:<%s}' % nameL).format(name) + rdict[name] + '\n'
        info.append(line)


    info = ''.join(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)



