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
        print("Usage: *.py count_sam.stat", file=sys.stderr)
        sys.exit()

    adict = {}
    with open(fn) as fp:
        for line in fp:
            name, read_count, char_count = line.strip().split()
            sp = cmn.lastName(name).split('_')[0]
            try:
                current = adict[sp]
            except KeyError:
                current = [0, 0]

            new = [current[0] + int(read_count), current[1] + int(char_count)]
            adict[sp] =new


    header = 'sp read_count char_count'.split()
    table = ['\t'.join(header)]
    for sp in adict:
        a, b = adict[sp]
        table.append('%s\t%s\t%s' % (sp, a, b))

    table.append('')
    print('\n'.join(table))




