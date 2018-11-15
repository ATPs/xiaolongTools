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
        wdir = sys.argv[1]
    except:
        print("Usage: *.py wdir", file=sys.stderr)
        sys.exit()

    total = 0
    unaligned = 0
    bestDir = cmn.txt_read('%s/best_mapping.txt' % wdir)
    fns = cmn.cmd2lines('ls %s/%s/*.sam' % (wdir, bestDir))
    sp = fns[0].split('/')[-3]
    for fn in fns:
        fp = open(fn)
        #with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                continue
            else:
                total += 1
                if line.strip().split()[2] == '*':
                    unaligned += 1
        fp.close()
    print(sp, bestDir, (total - unaligned), total)




