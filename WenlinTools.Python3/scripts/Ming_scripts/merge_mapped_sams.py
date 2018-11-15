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

def filter_and_write_lines(fp_dn, fn, header=False):
    with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                if header:
                    fp_dn.write(line)
            else:#alignment
                if line.strip().split()[2] != '*':
                    fp_dn.write(line)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def merge_sams(dn, fns):
    #dn = '%s.sam' % label

    print('merging files: %s into %s' % (str(fns), dn))

    if cmn.filexist(dn):
        cmn.run('rm ' + dn)

    fp_dn = open(dn,"a")

    filter_and_write_lines(fp_dn, fns[0], header=True)


    for fn in fns[1:]:
        filter_and_write_lines(fp_dn, fn)

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




