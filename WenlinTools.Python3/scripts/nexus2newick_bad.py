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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def remove_annotation(tree):
    new = []
    taken = True
    for char in tree:
        if char == '[':
            taken = False
        elif char == ']':
            taken = True
            continue

        if taken:
            new.append(char)
    return ''.join(new)

def replace_names(tree, mapdict):
    for name in mapdict:
        tree = tree.replace('(' + name + ':', '(' + mapdict[name] + ':')
    return tree


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    tFlag = False

    mapdict = {}
    trees = []
    for line in cmn.file2lines(fn):
        if 'Translate' in line:
            tFlag = True
            continue
        elif tFlag and ';' in line:
            tFlag = False

        if tFlag:
            key, name = line.strip(' ,').split()
            mapdict[key] = name


        if line.startswith('tree'):
            tree = '='.join(line.split('=')[1:]).strip()
            trees.append(tree)

    for tree in trees:
        tree = remove_annotation(tree)
        tree = replace_names(tree, mapdict)
        print(tree + '\r')

