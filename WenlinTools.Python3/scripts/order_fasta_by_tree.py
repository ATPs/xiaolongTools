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
        fn, ftree = sys.argv[1:3]
    except:
        print("Usage: *.py fasta ftree", file=sys.stderr)
        sys.exit()


    #ftree = '/project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/current_tree.newick'

    t = ete3.Tree(cmn.txt_read(ftree))

    order_list = []
    nameDict = {}
    for node in t:
        name = node.name.split('_')[0].lstrip("'")
        #name = node.name
        print(name)
        #if '_cp1' in node.name:
        order_list.append(name)
        nameDict[name] = node.name


    #read in fasta
    seqDict = read_fa(fn)
    mapDict = {key.split('_')[0]: key for key in seqDict}

    new = []
    taken_list = set([])
    missing_list = set([])
    for name in order_list:
        nodeName = nameDict[name]
        sp = name.split('_')[0]
        try:
            seqName = mapDict[sp]
            seq = seqDict[seqName]
            taken_list.add(name)
        except KeyError:
            #no such sequence
            missing_list.add(name)
            continue
        fasta = '>%s\n%s\n' % (nodeName, seq)
        new.append(fasta)

    dn = 're-ordered_seq.fa'
    cmn.write_file(''.join(new), dn)

    if len(missing_list) != 0:
        print('these following IDs in tree don\'t have sequences')
        print('\n'.join(missing_list))


    aset = set(seqDict.keys()) - taken_list
    if len(aset) != 0:
        print('the following IDs didn\'t show up in the provided tree')
        print('\n'.join(aset))
