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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()

    gap_chars = ['N', '-']
    Ngap = 0
    Ntotal = 0
    Ncount = 0
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line.strip()[1:]
                Ncount += 1
            else:
                seq = line.strip()
                for char in gap_chars:
                    Ngap += seq.count(char)
                Ntotal += len(seq)
    print(fn, Ncount, Ngap, Ntotal, float(Ngap) / Ntotal)



