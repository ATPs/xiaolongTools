import cmn
import sys


def aligned_reads(fn):
    #unaligned = 0
    #total = 0
    #previousID = ''
    aligned = {}
    unaligned = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                continue
            
            items = line.strip().split()
            #if items[0] == previousID:
            #    previousID = items[0]
            #    continue
            #total += 1
            ID = items[0]
            if items[2] == '*':
                try:
                    unaligned[ID] += 1
                except KeyError:
                    unaligned[ID] = 1
            else:
                try:
                    aligned[ID] += 1
                except KeyError:
                    aligned[ID] = 1
            
    #tell if the reads are paired
    if max(aligned.values()) > 1:#paired reads
        alignN = 0
        for ID in aligned:
            N = aligned[ID]
            if N > 2:
                alignN += 2
            else:
                alignN += N
        
        unalignN = sum(unaligned.values())

        total = alignN + unalignN
    else:
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

