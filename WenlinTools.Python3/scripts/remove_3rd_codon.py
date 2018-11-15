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
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    dn = cmn.lastName(fn).replace('.fa', '') + '_no3rd.fa'

    with open(fn) as fp, open(dn, 'w') as dp:
        for line in fp:
            if line[0] == '>':
                dp.write(line)
            else:
                seq = line.strip()
                for i, char in enumerate(seq):
                    if i % 3 == 2:
                        continue
                    dp.write(char)
                dp.write('\n')




