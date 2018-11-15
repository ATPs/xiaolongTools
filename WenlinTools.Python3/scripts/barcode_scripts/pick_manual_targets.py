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
        print("Usage: *.py *.report", file=sys.stderr)
        sys.exit()

    fcomp = 'compare.check'

    cmd = 'grep -v same compare.check'
    lines = cmn.cmd2lines(cmd)

    #exclude those with gap0 and nothing in report file
    badSp = set([])
    for line in cmn.file2lines(fn):
        if line[0] == '#':
            continue
        sp = line.split()[0].split('_')[0]
        badSp.add(sp)


    for line in lines:
        items = line.strip().split()
        sp = items[0]
        if sp not in badSp and items[2] == 'Gap0':
            if 'no' in items[-1]:
                line = 'goodReport:' + line
                print(line)
            continue
        print(line)


