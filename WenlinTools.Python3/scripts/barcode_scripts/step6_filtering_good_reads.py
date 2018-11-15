#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
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
        fn, fg = sys.argv[1:3]
    except:
        print("Usage: *.py 2nd_sam_aln.txt good_reads.txt", file=sys.stderr)
        sys.exit()

    goodIDs = set(cmn.getid(fg))

    dp = open('filtered_sam_aln.txt', 'w')
    dbad = open('bad_sam_aln.txt', 'w')
    with open(fn) as fp:
        for line in fp:
            Id = line.strip().split()[0]
            if Id in goodIDs:
                dp.write(line)
            else:
                dbad.write(line)
    dp.close()
    dbad.close()






