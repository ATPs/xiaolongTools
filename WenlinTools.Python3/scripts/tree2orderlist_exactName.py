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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        ftree = sys.argv[1]
    except:
        print("Usage: *.py ftree", file=sys.stderr)
        sys.exit()


    #ftree = '/project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/current_tree.newick'

    t = ete3.Tree(cmn.txt_read(ftree))
    # Calculate the midpoint node
    #R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    #t.set_outgroup(R)

    order_list = []
    nameDict = {}
    for node in t:
        name = node.name
        order_list.append(name)

    print('\n'.join(order_list))
