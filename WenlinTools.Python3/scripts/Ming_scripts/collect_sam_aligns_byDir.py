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

def count_sam_align(fns):
    total = 0
    unaligned = 0
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

    items = fn.split('/') 
    sp, ref = items[-3:-1]
    print(sp, ref, (total - unaligned), total)





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        wdir = sys.argv[1]
    except:
        print("Usage: *.py 3935", file=sys.stderr)
        sys.exit()

    #subdirs = cmn.cmd2lines('ls %s| grep -v txt$' % wdir)

    #for subdir in subdirs:
    #    fns = cmn.cmd2lines('ls %s/%s/*.sam' % (wdir, subdir))
    #    count_sam_align(fns)
    fn = '%s/mapping_stat.txt' % wdir
    if not cmn.filexist(fn):
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/count_sam_aligns_byDir.py %s' %  wdir
        cmn.run(cmd)
        sys.exit()

    sp = cmn.lastName(wdir)
    for line in cmn.file2lines(fn):
        print('%s\t%s' % (sp, line))



