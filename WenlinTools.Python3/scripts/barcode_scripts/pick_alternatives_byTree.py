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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_closest_node_name(node_list, good_names, target_node, t):
    order_nodes = sorted(node_list, key=lambda x: target_node.get_distance(x))
    for node in order_nodes:
        name = node.name
        genus = name2genus(name)
        if genus in good_names:
            return genus
    return None


def find_tree_node(t, ID):
    for node in t:
        if node.name.split('_')[0] == ID:
            return node
    return None


def name2genus(name):
    genus = name.replace('_', ' ').split()[1]
    return genus

if __name__=='__main__':
    #options=parse_options()
    try:
        fn, ftree = sys.argv[1:]
    except:
        print("Usage: *.py todoIDs ftree.renamed", file=sys.stderr)
        sys.exit()

    all_deflines = cmn.cmd2lines('grep ">" /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa|cut -d ">" -f 2')

    all_genus = set([each.split('_')[0] for each in all_deflines])

    t = ete3.Tree(ftree)
    # Calculate the midpoint node
    R = t.get_midpoint_outgroup()
    # and set it as tree outgroup
    t.set_outgroup(R)

    for line in cmn.file2lines(fn):
        ID = line.strip().split()[0]
        target_node = find_tree_node(t, ID)
        if target_node == None:
            print('%s is not in the tree' % ID)
            continue
        isFound = False
        node = target_node
        node = node.up
        while(not isFound and node != None):
            node_list = [each for each in node]
            names = [each.name for each in node_list]
            genus_list = set([name2genus(each) for each in names])
            good_names = genus_list & all_genus
            if len(good_names) > 0:
                isFound = True
                closest_name = find_closest_node_name(node_list, good_names, target_node, t)
            else:
                node = node.up
        print(ID, closest_name)



