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
def find_sp(line):
    sp = line[1:].split('_')[0]
    if ']' in sp:
        sp = sp.split(']')[-1]
    return sp


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    dn = 'labeled.fa'

    badIDs = set(cmn.file2lines('badIDs'))

    taken = set([])
    with open(fn) as fp, open(dn, 'w') as dp:
        for line in fp:
            if '(assembled)' not in line:
                dp.write(line)
            else:
                sp = find_sp(line)
                if sp in badIDs:
                    print('label %s as notSure' % sp)
                    line = '>[notSure]' + line[1:]
                    taken.add(sp)
                    dp.write(line)
                else:
                    dp.write(line)

    missing = badIDs - taken

    if len(missing) != 0:
        print('These badIDs are not labeled in the file:')
        print('\n'.join(missing))


