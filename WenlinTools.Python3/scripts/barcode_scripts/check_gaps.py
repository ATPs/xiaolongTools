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
        fn=sys.argv[1]
    except:
        print("Usage: *.py fa", file=sys.stderr)
        sys.exit()

    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4  in set([2,3]):
                continue

            if line[0] == '>':
                ID = line[1:]
            else:
                if 'N' in line or ('-' in line) or ('X' in line):
                    print(ID)
                    continue

                if any([char.islower() for char in line]):
                    print('lower', ID)





