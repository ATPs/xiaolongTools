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

#import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #fn = 'coding.fasta'
    fn = sys.argv[1]

    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line.strip() == '':
                continue
            if line[0] == '>':
                defline = line.strip()
            else:
                seq = line.strip()
                adict[defline] = seq

    length = len(seq)

    #informative positions are the positions that
    #at least group all the squence into 2 groups
    goodP = 0
    gapChars = set('X - N'.split())
    for i in range(length):
        groups = {}
        for name in adict:
            key = name.split('_')[0]
            char = adict[name][i]
            if char in gapChars:
                continue

            try:
                groups[char].add(key)
            except KeyError:
                groups[char] = set([key])
        group_count = 0
        for char in groups:
            aset = groups[char]
            if len(aset) > 1:
                group_count += 1

        if group_count >= 2:
            goodP += 1

    print('%s\t%s\t%s' % (fn, goodP, length))
