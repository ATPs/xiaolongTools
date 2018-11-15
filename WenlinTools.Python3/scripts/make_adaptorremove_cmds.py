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
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_adaptor_info():
    f1 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/adaptor/index_and_adaptor'
    bdict = {}
    with open(f1) as fp:
        for line in fp:
            items = line.strip().split()
            bdict[items[1]] = items[2]


    fdict = {}
    f2 = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/adaptor/lane_index.txt'
    with open(f2) as fp:
        for line in fp:
            name, index = line.strip().split()
            fdict[name] = bdict[index]
    return fdict



def group_reads(fns):
    adict = {}
    for fn in fns:
        fn = os.path.abspath(fn)
        label = '_'.join(cmn.lastName(fn).split('_')[:-1])
        try:
            adict[label].append(fn)
        except:
            adict[label] = [fn]
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py filelist", file=sys.stderr)
        sys.exit()


    fns = cmn.getid(fn)

    group_dict = group_reads(fns)

    adaptor_dict = read_adaptor_info()

    cmds = []
    for label in group_dict:
        #label = '_'.join(cmn.lastName(fn).split('_')[:2])
        fns = group_dict[label]
        if len(fns) != 2:
            print('number of pairs is incorrect for: %s' % str(fn))
        adaptor = adaptor_dict[label]
        cmd = '/home2/wli/local/AdapterRemoval-1.5.4/AdapterRemoval --file1 %s --file2 %s --basename %s ' % (fns[0], fns[1], label)
        cmd += '--trimns --pcr1 %s --pcr2 AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT ' % adaptor
        cmds.append(cmd)

    dn = 'ad_remove.cmds'
    cmn.write_lines(cmds, dn)



