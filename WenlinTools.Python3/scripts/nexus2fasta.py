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

    start = False
    new = []
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            if line == '':
                continue

            if line.lower().startswith('matrix'):
                start = True
                continue

            if start:
                try:
                    name, seq = line.split()
                except:
                    break
                new.append('>%s\n%s\n' % (name, seq))

    dn = cmn.lastName(fn) + '.fa'
    cmn.write_file(''.join(new), dn)




