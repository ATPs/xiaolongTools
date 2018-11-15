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
    rdict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        #sp = defline.split('_')[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, rdict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        ftree = sys.argv[2]
    except:
        print("Usage: *.py fa ftree", file=sys.stderr)
        sys.exit()


    #make the tree

    #ftree = '/project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/current_tree.newick'

    t = ete3.Tree(cmn.txt_read(ftree))

    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    order_list = []
    nameDict = {}
    for node in t:
        name = node.name.strip('\'"')
        #if '_cp1' in node.name:
        if True:
            order_list.append(name)
            #nameDict[name] = node.name.replace('_cp1', '').replace('/', '_')

    print('order_list', order_list)
    #read in fasta
    seqDict, rename_dict = read_fa(fn)

    new = []
    taken_list = set([])
    missing_list = set([])
    for name in order_list:
        #nodeName = rename_dict[name]
        nodeName = name
        try:
            seq = seqDict[name]
            taken_list.add(name)
        except KeyError:
            #no such sequence
            missing_list.add(name)
            continue
        fasta = '>%s\n%s\n' % (nodeName, seq)
        new.append(fasta)

    if len(missing_list) != 0:
        print('these following IDs in tree don\'t have sequences')
        print('\n'.join(missing_list))


    aset = set(seqDict.keys()) - taken_list
    if len(aset) != 0:
        print('the following IDs didn\'t show up in the provided tree')
        print('\n'.join(aset))
        for name in aset:
            fasta = '>%s\n%s\n' % (name, seqDict[name])
            new.append(fasta)

    dn = 're-ordered_seq.fa'
    cmn.write_file(''.join(new), dn)

