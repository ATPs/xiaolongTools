#!/home2/wli/anaconda/bin/python

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
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def aligned_reads(fn):
    #unaligned = 0
    #total = 0
    #previousID = ''
    aligned = set([])
    unaligned = set([])
    samfile = pysam.AlignmentFile(fn)
    for record in samfile:
        ID = record.query_name
        if record.is_read1:
            ID += '/1'
        elif record.is_read2:
            ID += '/2'
        
        if record.is_unmapped:
            unaligned.add(ID)
        else:
            aligned.add(ID)
    
    alignN = len(aligned)
    total = alignN + len(unaligned)

    return alignN, total


def count_sam_align(fns):
    totalN = 0
    alignN = 0
    for fn in fns:
        alnN, tN = aligned_reads(fn)
        totalN += tN
        alignN += alnN

    items = fn.split('/') 
    sp, ref = items[-3:-1]
    print(sp, ref, alignN, totalN)





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

    subdirs = cmn.cmd2lines('ls %s| grep -v txt$' % wdir)

    for subdir in subdirs:
        fns = cmn.cmd2lines('ls %s/%s/*.sam' % (wdir, subdir))
        count_sam_align(fns)



