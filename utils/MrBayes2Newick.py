
import re
from ete3 import Tree

def readMrBayesTree(file_in, file_out=None):
    '''
    given a output of mrbayes, return a str of tree in newick format, with prob(percent) as bootstrap values
    if file_out is provided, then save the file to file_out
    '''
    txt = open(file_in).read()
    txt = txt.replace('\n','')
    tree_pattern = re.compile('begin trees;.*?end;', re.M|re.I)
    tree = tree_pattern.findall(txt)
    tree = tree[0].strip(';')
    elements = tree.split(';')
    txt_species = elements[1]
    txt_species = txt_species.split()
    dc_species = {}
    for n in range(1,len(txt_species),2):
        dc_species[txt_species[n]] = txt_species[n+1].strip(',')
    
    txt_tree = elements[2]
    txt_tree = txt_tree.split(' ')[-1]
    ls_tree_info = re.split('\[|\]',txt_tree)
    ls_tree_info_keep = []
    for n,e in enumerate(ls_tree_info):
        if e.startswith('&'):
            if 'prob(percent)=' in e:
                e = e.split('prob(percent)="')[1].split('"')[0]
                if ls_tree_info[n-1].endswith(')'):
                    ls_tree_info_keep.append(e)
        else:
            ls_tree_info_keep.append(e)
    result_tree = ''.join(ls_tree_info_keep) +';'
    tree = Tree(result_tree)
    for leaf in tree.iter_leaves():
        leaf.name = dc_species[leaf.name]
    
    tree_txt = tree.write()
    if file_out is not None:
        tree.write(outfile=file_out) 
    
    return tree_txt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='input a MrBayes result, output a newick file. Only keep the probability at the branch sites as bootstrap values')
    parser.add_argument('-i','--input', help = 'input file of MrBayes result', required=True)
    parser.add_argument('-o','--output',help = 'location of the where the output file stored',default=None)
    f = parser.parse_args()
    readMrBayesTree(file_in=f.input, file_out=f.output)