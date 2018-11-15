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
        print("Usage: *.py samfile", file=sys.stderr)
        sys.exit()

    total = 0
    aligned = 0
    new = []
    with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                new.append(line)
                continue
            else:
                total += 1
                if line.strip().split()[2] != '*':
                    aligned += 1
                    new.append(line)

    print(cmn.lastName(fn), aligned, total)
    dn = fn.replace('.sam', '_mapped.sam')
    cmn.write_file(''.join(new), dn)

    print('result is in %s' % dn)




