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
        Ncores = int(sys.argv[2])
    except:
        print("Usage: *.py fa Ncores", file=sys.stderr)
        sys.exit()

    #if the nodes are less than 4 taxa, produce a random tree
    cmd = "grep '>' %s" % (fn)
    lines = [each[1:].strip() for each in cmn.cmd2lines(cmd)
            if each.strip() != '']

    N = len(lines)
    if N < 4:
        print('Warning: fastme can not make tree of less than 4 taxa')
        print('Warning: so I make a fake tree...')
        dn = '%s.phylip.fastme.tre' % cmn.lastName(fn)
        if N == 1:
            info = '(%s);\n' % lines[0]
        if N == 2:
            a, b = lines
            info = '(%s,%s);\n' % (a, b)
        elif N == 3:
            a, b, c = lines
            info = '((%s,%s),%s);\n' % (a, b, c)
        cmn.write_file(info, dn)
        sys.exit()

    label = cmn.lastName(fn)
    cmd = 'rm RAxML_*.%s;' % label
    cmd += '/home2/wli/local/RAxML/raxmlHPC-PTHREADS-SSE3 -m GTRGAMMA -p 7112 -T %s -s %s -n %s' % (Ncores, label, label)
    cmn.run(cmd)




