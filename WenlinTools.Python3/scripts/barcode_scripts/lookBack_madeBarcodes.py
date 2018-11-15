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
        sp = defline.split('_')[0]
        seq = ''.join(lines[1:])
        adict[sp] = seq
        rdict[sp] = defline
    return adict, rdict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_madePaired(fn, label=None):
    adict = {}
    chunk = []
    found_pair = False
    fp = open(fn)
    for line in fp:
        if '(assembled)' in line:
            #check if it is a new line
            if found_pair and len(chunk) != 0:
                #report the last
                adict[sp] = ''.join(chunk)
                chunk = []
                found_pair = False

            sp = line[1:].split('_')[0].split(']')[-1]

            if label != None:
                line = '>%s%s' % (label, line[1:])

        elif line[0] == '>':
            #has def line and no 'assembled'
            found_pair = True

        chunk.append(line)
    fp.close()

    adict[sp] = ''.join(chunk)
    return adict



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        ftree = sys.argv[2]
    except:
        print("Usage: *.py paired.fasta nucleus.tre", file=sys.stderr)
        sys.exit()


    fpair = '/project/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/paired_verified_barcodes.fa'
    madeDict = parse_madePaired(fpair, label='[old]')
    #ftree = '/project/biophysics/Nick_lab/wli/sequencing/Eudamine/BEAST_timing/current_tree.newick'

    newDict = parse_madePaired(fn, label='[new]')

    allDict = dict(madeDict)
    allDict.update(newDict)

    t = ete3.Tree(cmn.txt_read(ftree))

    # Calculate the midpoint node
    #R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    #t.set_outgroup(R)

    order_list = []
    nameDict = {}
    for node in t:
        sp = node.name.split('_')[0].lstrip("'")
        if '_cp1' in node.name:
            order_list.append(sp)
            #nameDict[sp] = node.name.replace('_cp1', '').replace('/', '_')
    #print order_list
    #print allDict.keys()

    #read in fasta
    #seqDict, rename_dict = read_fa(fn)

    new = []
    taken_list = set([])
    missing_list = set([])
    for sp in order_list:
        #nodeName = nameDict[name]
        try:
            fasta = allDict[sp]
            taken_list.add(sp)
            if sp in madeDict:
                print('addBack: %s' % sp)

        except KeyError:
            missing_list.add(sp)
            continue
        new.append(fasta)

    if len(missing_list) != 0:
        print('these following IDs in tree don\'t have sequences')
        print('\n'.join(missing_list))


    aset = set(newDict.keys()) - taken_list
    if len(aset) != 0:
        print('the following IDs didn\'t show up in the provided tree')
        print('\n'.join(aset))
        for sp in aset:
            fasta = newDict[sp]
            new.append(fasta)

    dn = 'addBack_seq.fa'
    cmn.write_file(''.join(new), dn)

