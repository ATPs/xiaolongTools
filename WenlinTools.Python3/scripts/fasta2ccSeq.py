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
gapChars = set(['N', '-', 'X'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        #defline, seq = each.strip().split()
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:]).replace('N', '-')
        adict[defline] = seq
    return adict, len(seq)


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

    seqDict, length = read_fa(fn)

    seqs = list(seqDict.values())
    new = []
    for i in range(length):
        chars = [seq[i] for seq in seqs]
        chars = [char for char in chars
                if char not in gapChars]
        if len(chars) == 0:
            new.append('-')
        else:
            count_dict = Counter(chars)
            char = max(count_dict, key=lambda x: count_dict[x])
            if count_dict[char] >= (0.4 * len(chars)):
                new.append(char)
            else:
                new.append('N')

    print('>ccSeq\n%s\n' % ''.join(new))


