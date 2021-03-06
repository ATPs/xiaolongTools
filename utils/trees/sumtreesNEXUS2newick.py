# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 12:07:49 2018

@author: ATPs
"""
import dendropy
from ete3 import Tree

def treeNexus2newick(filein,fileout=None):
    '''
    filein is a filename of a tree in NEXUS format, eg. generated by summtrees
    write a tree to fileout in newick format
    if fileout is None, fileout = filein +'.nwk'
    '''
    if fileout is None:
        fileout = filein + '.nwk'
    

    mle = dendropy.Tree.get(path=filein, schema="nexus")
    tree_txt = mle.as_string(schema='newick')
    tree = Tree(tree_txt.split()[1].strip())
    tree.resolve_polytomy(recursive=True)
    tree.write(outfile=fileout)

description = '''filein is a filename of a tree in NEXUS format, eg. generated by summtrees
write a tree to fileout in newick format
if fileout is None, fileout = filein +'.nwk'
'''
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='convert a tree from NEXUS format to newick')
    parser.add_argument('-i','--input', help = 'input newick tree file', required=True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored',default=None)
    f = parser.parse_args()
    treeNexus2newick(filein=f.input, fileout=f.output)