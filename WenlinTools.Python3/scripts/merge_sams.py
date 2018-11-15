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

def merge_sams(dn, fns):
    #dn = '%s.sam' % label

    print('merging files: %s into %s' % (str(fns), dn))

    if cmn.filexist(dn):
        cmn.run('rm ' + dn)

    cmn.run('cp %s %s' % (fns[0], dn))

    fp_dn = open(dn,"a")
    for fn in fns[1:]:
        fp = open(fn)
        for line in fp:
            #if line[0] != "@" and line[0] != "[" and line.split()[2] != "*":
            if line[0] != "@":
                fp_dn.write(line)
        fp.close()

    fp_dn.close()
    return dn



if __name__=='__main__':
    #options=parse_options()
    try:
        dn = sys.argv[1]
        fns=sys.argv[2:]
    except:
        print("Usage: *.py dn sams_to_merge", file=sys.stderr)
        sys.exit()

    merge_sams(dn, fns)




