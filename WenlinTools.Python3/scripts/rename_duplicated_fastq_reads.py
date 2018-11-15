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
        print("Usage: *.py fq", file=sys.stderr)
        sys.exit()

    check_set = set([])

    new = []
    rename_label = 0
    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 0:
                items = line.strip().split()
                ID = items[0]
                if ID in check_set:
                    items[0] += '/V%s' % rename_label
                    rename_label += 1
                    line = ' '.join(items) + '\n'
                else:
                    check_set.add(ID)
            new.append(line)

    dn = fn + '.renamed'
    cmn.write_file(''.join(new), dn)





