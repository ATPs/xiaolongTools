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
import random


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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
                label = line.strip()
            else:
                seq = line.strip()
                adict[label] = '%s\n%s\n' % (label, seq)


    times = 10
    keys = list(adict.keys())
    cmn.mkdir('shuffle_genome')
    for each in range(times):
        random.shuffle(keys)
        new = [adict[key] for key in keys]
        dn = 'shuffle_genome/%s_shuffle%s' % (cmn.lastName(fn), each)
        cmn.write_file(''.join(new), dn)



