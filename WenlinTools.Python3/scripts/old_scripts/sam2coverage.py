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
from collections import Counter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rdict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        '-': 'N',
        '*': 'N'
        }

def reverse_seq(seq):
    return ''.join([rdict[char] for char in seq[::-1]])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.sam", file=sys.stderr)
        sys.exit()

    record_dict = {}
    with open(fn) as fp:
        for line in fp:
            if line.startswith('@'):
                continue
            items = line.strip().split('\t')
            position = int(items[3])
            seq = items[9]
            qual = items[10]
            scaffold = items[2]
            flag = items[1]
            #bits = bin(int(flag)).replace('0b', '')[::-1]
            #if ((len(bits) >= 5 and bits[4] == '1')
            #        or (len(bits) >= 6 and bits[5] == '1')):
            #    seq = reverse_seq(seq)
            #if len(bits) >= 3 and bits[2] == '1':
                #print 'skip quality flag %s' % flag
            #    continue
            if scaffold not in record_dict:
                record_dict[scaffold] = {}
            for i, char in enumerate(seq):
                ii = position + i
                try:
                    record_dict[scaffold][ii].append(char)
                except KeyError:
                    record_dict[scaffold][ii] = [char]

    dn = fn + '.coverage'
    new = []
    for scaffold in record_dict:
        position_dict = record_dict[scaffold]
        for position in position_dict:
            counts = position_dict[position]
            freq_dict = Counter(counts)
            stat = '|'.join(['%s:%s' % (char, freq_dict[char]) for char in freq_dict])
            new.append('%s\t%s\t%s' % (scaffold, position, stat))

    cmn.write_lines(new, dn)


