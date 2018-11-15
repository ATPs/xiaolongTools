#!/home2/wli/anaconda/bin/python

import cmn
import sys
import pysam


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


def old_aligned_reads(fn):
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


wdir = sys.argv[1]

subdirs = cmn.cmd2lines('ls %s| grep -v txt$' % wdir)

if len(subdirs) == 1:
    print('only one directory, no need to tell who is best')
    cmn.write_file(subdirs[0], '%s/best_mapping.txt' % wdir)
    sys.exit()


else:
    sizeDict = {}
    for each in subdirs:
        fsams = cmn.cmd2lines('ls %s/%s/*.sam' % (wdir, each))
        total = 0
        mapped = 0
        halfmap = 0
        TpN = 0
        TptN = 0
        for fsam in fsams:
            if not cmn.filexist(fsam):
                continue
            mappedN, halfN, totalN, pN, ptN = aligned_reads(fsam)
            total += totalN
            mapped += mappedN
            halfmap += halfN
            TpN += pN
            TptN += ptN

        sizeDict[each] = (mapped, total, halfmap, float(TpN)/TptN)

    dn = '%s/mapping_stat.txt' % wdir
    info = ['%s\t%s\n' % (name, '\t'.join(map(str, sizeDict[name])))
            for name in sizeDict]
    cmn.write_file(''.join(info), dn)

    best = max(list(sizeDict.keys()), key = lambda x: sizeDict[x][2])
    print('best mapping', best)

    cmn.write_file(best, '%s/best_mapping.txt' % wdir)

