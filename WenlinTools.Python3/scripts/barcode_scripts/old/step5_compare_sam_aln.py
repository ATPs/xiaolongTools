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

def parse_bait_file(fbait):
    new = {}
    target = {}
    for line in cmn.file2lines(fbait):
        ID, sp, seq = line.strip().split()
        #name = '%s_%s' % (ID, sp)
        name = sp
        if len(target) == 0:
            new[name] = seq

        target[name] = seq
    return new, target

def get_top_and_target_barcodes(fn, topN):
    seqDict = read_fa('/work/00412/mtang/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa')

    count_dict = {}
    for line in cmn.file2lines(fn):
        count, name = line.strip().split()
        count_dict[name] = count

    names = sorted(list(count_dict.keys()), key=lambda x: count_dict[x], reverse=True)

    new = {}
    target = {}
    for i in range(topN):
        name = names[i]
        if i == 0:
            target[name] = seqDict[name]
        new[name] = seqDict[name]
    return new, target

def read_aln(fn):
    seqs = {}
    for line in cmn.file2lines(fn):
        name, seq = line.strip().split()
        seqs[name] = seq

    return seqs

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict

def name2sp(name):
    sp = '('.join(name.split('|')[-1].split('(')[:-1]).strip('_')
    return sp

def detect_mismatch(seq, read_alg):
    #consider every letter, even those not aligned by bwa
    M = 0
    length = len(seq)
    for i, char in enumerate(read_alg):
        char = char.upper()
        if char == '-' or char == 'N':
            continue
        if i >= length:
            continue

        char2 = seq[i].upper()
        if char2 == '-' or char2 == 'N':
            continue

        if char != char2:
            M += 1
    return M

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def compute_identity(a, b):
    print(a)
    print(b)
    gapChars = set(['-', 'N'])
    check_list = [a[i] == b[i] for i in range(len(a))
            if a[i] not in gapChars and b[i] not in gapChars]
    identity = float(sum(check_list)) / len(check_list)
    return identity


if __name__=='__main__':
    #options=parse_options()
    try:
        fn1, fn2, fbait = sys.argv[1:4]
    except:
        print("Usage: *.py 1stSam 2ndSam sampleInfo.bait", file=sys.stderr)
        sys.exit()

    fn = '/work/00412/mtang/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #1. compare ID to see overlapping and uniq ID
    barcodeDict = read_fa(fn)
    #refDict, targetDict = get_top_and_target_barcodes(fcount, 4)
    refDict, targetDict = parse_bait_file(fbait)

    identity_cut = min([compute_identity(refDict[name1], refDict[name2]) for name1 in refDict for name2 in refDict])
    identity_cut -= 0.1
    print('identity_cut', identity_cut)

    seqDict1 = read_aln(fn1)
    seqDict2 = read_aln(fn2)
    #refSp = name2sp(next(iter(seqDict2)))

    #ID1 is the mapping with expected bait
    #ID2 is the IDs to exclude contamination
    #So good ID should be the ones in ID1 but not in ID2
    ID1mapping = {name.split('|')[0].split('(')[0] : name for name in seqDict1}
    ID2mapping = {name.split('|')[0].split('(')[0] : name for name in seqDict2}

    ID1 = set(ID1mapping.keys())
    ID2 = set(ID2mapping.keys())
    #refID1 = set([ID for ID in ID1mapping
    #        if name2sp(ID1mapping[ID]) in refDict])
    #for ID in refID1:
    #    if ID not in ID2mapping:
    #        ID2mapping[ID] = ID1mapping[ID]
    #things only in ID1 and not in ID2 is really bad
    #bad_IDs = ID1 - ID2 - refID1
    #good_IDs = (ID2 - ID1) | refID1
    #bad_IDs = ID1 - ID2
    bad_IDs = set([])
    good_IDs = ID1 - ID2

    overlapIDs = ID1 & ID2

    #print 'bad IDs', len(bad_IDs)
    print('good IDs', len(good_IDs))
    print('overlapping IDs', len(overlapIDs))
    #print bad_IDs
    print(good_IDs)

    #refSeq = barcodeDict[refSp]
    refSp = list(targetDict.keys())[0]
    refSeq = targetDict[refSp]

    good_reads = []
    bad_reads = []
    for ID in overlapIDs:
        name1 = ID1mapping[ID]
        name2 = ID2mapping[ID]
        sp2 = name2sp(name2)
        if name2 in refDict:
            #the names are among the baits
            good_reads.append(name1)
            continue
        #1 is baited and 2 is bad
        #sp1 = name2sp(name1)
        print(ID, name1, name2)
        seq2 = barcodeDict[sp2]
        print(sp2, seq2, refSeq)
        identity = compute_identity(seq2, refSeq)

        #M1 is mismatch between ref
        aln1 = seqDict1[name1]
        misM1 = detect_mismatch(refSeq, aln1)

        #M2 is mismatch between bad IDs
        aln2 = seqDict2[name2]
        misM2 = detect_mismatch(seq2, aln2)

        print('comp: %s_|_%s = (%s:%s) %s' % (name1, name2, misM1, misM2, identity))

        if misM1 == misM2:
            #has equal mismatch, let's put it as True
            #this is a little dangous but we need more reads for those poor samples
            good_reads.append(name1)
        elif (misM2 + 1) >= misM1:
            #bad one has more mismatch, good one is good!
            #good one can be 1 bp more than the bad one
            good_reads.append(name1)
        else:
            if identity >= identity_cut:
                good_reads.append(name1)
            else:
                bad_reads.append((name2, aln1))

    print('further classify overlapping reads into:')
    print('%s good reads' % len(good_reads))
    print('%s bad reads' % len(bad_reads))
        #sp2 = name2sp(name2)

    #add back the previous IDs
    good_reads.append('#' * 100)
    for ID in good_IDs:
        name = ID1mapping[ID]
        good_reads.append(name)

    cmn.write_lines(good_reads, 'good_reads.txt')

    #bad_reads.append('#' * 100)
    for ID in bad_IDs:
        name = ID2mapping[ID]
        bad_reads.append((name, seqDict2[name]))
    bad_alignments = ['%s    %s\n' % (each[0], each[1]) for each in bad_reads]
    cmn.write_file(''.join(bad_alignments), 'bad_reads_alignment.txt')

