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
        for fsam in fsams:
            mappedN, totalN = aligned_reads(fsam)
            total += totalN
            mapped += mappedN
        sizeDict[each] = (mapped, total)

    dn = '%s/mapping_stat.txt' % wdir
    info = ['%s\t%s\t%s\n' % (name, sizeDict[name][0], sizeDict[name][1])
            for name in sizeDict]
    cmn.write_file(''.join(info), dn)

    best = max(list(sizeDict.keys()), key = lambda x: sizeDict[x][0])
    print('best mapping', best)

    cmn.write_file(best, '%s/best_mapping.txt' % wdir)

