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
from collections import Counter


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
            if line[0] == '>':
                defline = line.strip()
            else:
                seq = line.strip()
                adict[defline] = seq

    sequences = list(adict.values())
    length = len(sequences)

    final = {}
    for each in zip(*sequences):
        count_dict = Counter(each)
        maxCount = max(count_dict.values())
        try:
            gapF = count_dict['-']
        except:
            gapF = 0
        fraction = maxCount * 10 / (length - gapF)

        try:
            final[fraction] += 1
        except KeyError:
            final[fraction] = 1

    stat = 0
    for name in final:
        key = float(name) / 10
        value = final[name]
        print(key, value)
        if name != 10:
            stat += value

    #print 'total positions: %s' % len(sequences[0])
    #print 'informative positions: %s (%s)' % (stat, float(stat)/len(sequences[0]))
    print('%s\t%s\t%s' % (fn, len(sequences[0]), stat))
