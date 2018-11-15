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
import random

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def chars2line(pos, chars):
    pos = str(pos)
    line = [pos, pos]
    for i in range(0, len(chars), 2):
        a, b = chars[i:i+2]
        if a == b:
            line.append(a)
        else:
            line.append('%s,%s' % (a, b))
    return '\t'.join(line)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].split('_')[0]
            else:
                seq = line.strip()
                try:
                    adict[name].append(seq)
                except KeyError:
                    adict[name] = [seq]

    names = list(adict.keys())
    char_label = {}
    current_count = 1
    missing_data = set('N -'.split())
    #header = ['fake', 'position'] + names
    #header = '\t'.join(header)

    seqs = []
    #hetero = [[header], [header], [header], [header]]
    for name in names:
        seqs += adict[name]

    result = [['%s\t%s' % (name, i+1)] for name in names for i in range(2)]
    gapcut = 0.99
    #enrich more informative positions
    Ncut = 4 #each nucleatide appear should have larger count than this
    for i in range(len(seqs[0])):
        chars = [seq[i] for seq in seqs]
        check_chars = [char for char in chars if (char != '-') and (char != 'N')]
        if len(check_chars) < (gapcut * len(chars)):
            continue
        N = len(set(check_chars))
        count_dict = Counter(check_chars)
        if any([countN < Ncut for countN in list(count_dict.values())]):
            continue

        if N == 0 or N == 1:
            #skip the all gapped positions
            #also skip same character lines
            continue
        else:
            for i, char in enumerate(chars):
                if char in missing_data:
                    result[i].append(-9)
                    continue

                try:
                    code = char_label[char]
                except KeyError:
                    code = current_count
                    char_label[char] = current_count
                    current_count += 1

                result[i].append(code)
    #random sample
    p = 1000
    indexes = list(range(len(result[0])-1))
    goodI = random.sample(indexes, p)
    print(len(goodI))
    goodI.append(-1)

    goodI.sort()

    dn = cmn.lastName(fn) + 'STRUCTUREinput.txt'
    new = ['\t'.join(map(str, [line[index+1] for index in goodI])) for line in result]
    new.append('')
    cmn.write_lines(new, dn)




