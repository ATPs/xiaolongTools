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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, Range = sys.argv[1:3]
        i, j = list(map(int, Range.split('-')))
    except:
        print("Usage: *.py aln 0-10000", file=sys.stderr)
        sys.exit()

    new = []
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                new.append(line)
            else:
                seq = line.strip()
                new.append(seq[i:j] + '\n')

    dn = '%s_%s.fa' % (cmn.lastName(fn).replace('.fa', ''), Range)
    cmn.write_file(''.join(new), dn)



