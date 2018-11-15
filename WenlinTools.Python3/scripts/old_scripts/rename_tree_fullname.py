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

def get_names():
    adict = {}
    fns = cmn.cmd2lines('ls -tr /project/biophysics/Nick_lab/wli/sequencing/scripts/data/*.sampleData')
    for fn in fns:
        for line in cmn.file2lines(fn):
            line = line.strip()
            items = line.split()
            sp = items[0].split('-')[-1]
            line = line.replace(items[0], sp).replace('-', '_').replace('(','').replace(')', '')
            adict[sp] = '_'.join(line.split())
    return adict


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

    f_table = '/project/biophysics/Nick_lab/wli/sequencing/scripts/name_table'

    #nameDict = {}
    #for line in cmn.file2lines(f_table):
    #    items = line.strip().split()
    #    if len(items) == 0:
    #        continue
    #    label = items[0]
    #    name = '_'.join(items[1:])
    #    nameDict[label] = name.replace('-', '_')

    #print nameDict.keys()
    nameDict = get_names()

    t = ete3.Tree(cmn.txt_read(fn).replace('[&U]', ''))

    appear = {}
    for node in t:
        name = node.name
        sp = name.split('_')[0].split('.')[0]
        if sp not in appear:
            appear[sp] = 1
        else:
            appear[sp] += 1

        new_name = '%s_cp%s' % (nameDict[sp], appear[sp])
        node.name = new_name

    info = t.write()
    print(info)
    dn = fn + '.renamed'
    cmn.write_file(info, dn)



