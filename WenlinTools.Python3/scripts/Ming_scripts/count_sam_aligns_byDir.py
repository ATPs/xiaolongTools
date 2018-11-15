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
def map_fraction(record):
    aligned = record.get_aligned_pairs()
    N = 0
    alignN = 0
    for i, j in aligned:
        N += 1
        if i != None and j != None:
            alignN += 1
    
    fraction = float(alignN) / N
    return fraction


def aligned_reads(fn):
    #unaligned = 0
    #total = 0
    #previousID = ''
    aligned = set([])
    unaligned = set([])
    samfile = pysam.AlignmentFile(fn)
    pN, tN = 0, 0 #use to count stat by positions
    halfmap = set([]) # scaffolds aligned more than half positions
    for record in samfile:
        ID = record.query_name
        if record.is_read1:
            ID += '/1'
        elif record.is_read2:
            ID += '/2'
        
        if record.is_unmapped:
            unaligned.add(ID)
        else:#this is aligned
            aligned.add(ID)

            #check map fraction
            if map_fraction(record) >= 0.5:
                halfmap.add(ID)
        
        if not record.is_secondary:
            #count positions
            aligned_pairs = record.get_aligned_pairs()
            #i is read, j is ref
            for i, j in aligned_pairs:
                if i == None:
                    continue
                
                tN += 1
                if j != None:
                    pN += 1

    alignN = len(aligned)
    total = alignN + len(unaligned)
    half_alignN = len(halfmap)
    #pPercent = float(pN) / tN

    return alignN, half_alignN, total, pN, tN


def count_sam_align(fns):
    totalN = 0
    alignN = 0
    half_alignN = 0#more than half aligned
    total_pN = 0 #mapped positions
    total_ptN = 0 # total positions
    for fn in fns:
        if not cmn.filexist(fn):
            continue
        
        #pN and ptN are the counts by positions
        alnN, halfN, tN, pN, ptN = aligned_reads(fn)
        totalN += tN
        alignN += alnN
        half_alignN += halfN
        total_pN += pN
        total_ptN += ptN
    
    pPercent = float(total_pN) / total_ptN
    items = fn.split('/') 
    sp, ref = items[-3:-1]
    print(sp, ref, alignN, totalN, half_alignN, pPercent)




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



