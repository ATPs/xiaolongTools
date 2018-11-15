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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        seq = seq.replace('N', '-')
        adict[defline] = seq
    return adict


def print4order(line):
    items = line.replace('?','').replace('__', '_').strip().split('_')[1:]
    try:
        items.remove('F')
    except:
        pass
    try:
        items.remove('M')
    except:
        pass
    print(items)
    return items
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py fasta", file=sys.stderr)
        sys.exit()

    seqDict = read_fa(fn)

    names = list(seqDict.keys())

    sorted_names = sorted(names, key=lambda x: print4order(x))

    dn = 'sorted.fa'
    dp = open(dn, 'w')
    for name in sorted_names:
        fasta = '>%s\n%s\n' % (name, seqDict[name])
        dp.write(fasta)

    dp.close()
