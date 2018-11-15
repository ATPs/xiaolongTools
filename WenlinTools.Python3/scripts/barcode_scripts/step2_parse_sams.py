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
import pysam
import math

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bad_alignment = []

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def consensus_alignment(aligns):
    #actually if the reads didn't show conflict, we will take it
    ccAlign = {}
    #isGood = True
    indexes = set([])
    for each in aligns:
        print(each)

    for each in aligns:
        for i, j in each:
            #i is read index
            #j is ref
            indexes.add(i)
            if j == None:
                continue
            if i == None:
                continue

            try:
                ccAlign[i].append(j)
            except:
                ccAlign[i] = [j]
    indexes = list(indexes)
    indexes.sort()

    aligned = []
    alignedN = 0
    conflictN = 0
    consCut = 0.8
    for i in indexes:
        try:
            jList = ccAlign[i]
        except:
            jList = None

        if jList == None:
            #it is ok if no aligned
            aligned.append((i, None))
        else:#not None
            checkA = jList[0]
            #cutoff = math.ceil(consCut * len(jList))
            cutoff =consCut * len(jList)
            if jList.count(checkA) >= cutoff:
                aligned.append((i, checkA))
                alignedN += 1
            else:
                #it is not ok to show conflicts
                aligned.append((i, None))
                conflictN += 1
    #print 'aligned', alignedN, conflictN, aligned
    if alignedN == 0:
        return None
    elif float(conflictN) / alignedN <= 0.25:
        return aligned
    else:
        return None


def old_consensus_alignment(aligns):
    #actually if the reads didn't show conflict, we will take it
    ccAlign = {}
    isGood = True
    for each in aligns:
        for i, j in each:
            #i is read index
            #j is ref
            try:
                currentJ = ccAlign[i]
                if currentJ == None:
                    ccAlign[i] = j
                elif j == None:
                    continue
                elif currentJ != j:
                    print('diff', currentJ, j)
                    isGood = False
                    return None
            except:
                ccAlign[i] = j
        if not isGood:
            return None
    #if reach here, then the read is good
    #transform the align back
    aligned = []
    keys = list(ccAlign.keys())
    keys.sort()
    for i in keys:
        aligned.append((i, ccAlign[i]))
    return aligned


def number4sorting(seq):
    i = 0
    while (i < len(seq)):
        if not seq[i].isupper():
            i += 1
        else:#is upper now
            break
    j = i
    while (j < len(seq)):
        if seq[i].isupper():
            j += 1
        else:#has lower now
            break
    return (i, j)



def old_seqblock2alignment(aligned, seq):
    if True:
        subjct_range = [i[1] for i in aligned if i[1]!=None]
        left = min(subjct_range)
        #right = max(subjct_range)

        #change unmapped letter into lower case
        seq = list(seq)
        print(seq)
        print(aligned)
        for i, j in aligned:
            if i == None:
                continue
            if j == None:
                seq[i] = seq[i].lower()

        seq = ''.join(seq)


        if aligned[0][1] == None:
            #need to fill in beginning
            shifts = []
            for i, each in enumerate(aligned):
                if each[1] == None:
                    shifts.insert(0, -i - 1)
                else:
                    break

            print('before', left, aligned)
            for i, shift in enumerate(shifts):
                aligned[i] = (aligned[i][0], left + shift)
            left -= len(shifts)
            print('after', left, aligned)

            if left < 0:#shift too much
                for i in range(0 - left):
                    aligned[i] = (aligned[i][0], None)

                left = 0
                print('fix', left, aligned)


        #make the alignment output
        aln = []
        for _ in range(left):
            aln.append('-')

        isPass = False
        for i, j in aligned:
            if j == None:
                if isPass:
                    aln.append(seq[i])
                continue
            if i == None:
                aln.append('-')
                continue
            isPass = True
            aln.append(seq[i])

        aln = ''.join(aln)

        return aln


def group_by_spnames(name):
    key = '|'.join(name.split('|')[:-1])
    return key


def alnDict2output(aln_dict, dn, order='sorting'):
    info = []
    if len(aln_dict) == 0:
        cmn.run('touch %s' % dn)
        return None
    #maxLength = max([len(each) for each in aln_dict.keys()])
    maxLength = 0
    maxNameLength = max([len(each) for each in aln_dict])
    nameformat = '{:<%s}' % maxNameLength

    names = list(aln_dict.keys())
    if order == 'sorting':
        names = sorted(names, key=lambda x: number4sorting(aln_dict[x]))
    elif order == 'grouping':
        #this is used to output inconsistent group
        #rank by grouping of species IDs
        names = sorted(names, key=lambda x: group_by_spnames(x))
    else:
        names.sort()

    for i, name in enumerate(names):
        #name = 'readgroup%s' % i
        aln = aln_dict[name]
        name = nameformat.format(name)

        toAdd = maxLength - len(aln)
        if toAdd > 0:
            aln += '-' * toAdd
        info.append('%s    %s\n' % (name, ''.join(aln)))
    cmn.write_file(''.join(info), dn)


def record2name(record, showRef=False):
    name = record.query_name
    if record.is_read1:
        name += '_1'
    elif record.is_read2:
        name += '_2'
    if showRef:
        name = '%s|%s' % (name, record.reference_name)
    return name


def collapse_same_alignment(aln_dict):
    hashDict = {}
    for name in aln_dict:
        seq = aln_dict[name]
        try:
            hashDict[seq].append(name)
        except KeyError:
            hashDict[seq] = [name]

    new = {}
    for seq in hashDict:
        names = hashDict[seq]
        if len(names) == 1:
            new[names[0]] = seq
        else:
            newName = '[%sr]%s' % (len(names), '|'.join(names))
            new[newName] = seq
    return new

def seqblock2alignment(aligned, name, seq):
    #tranform the alinged into a dict for better reference
    align_dict = {}
    for i, j in aligned:
        if j != None:
            align_dict[j] = i

    #subjct_range = [i[1] for i in aligned if i[1]!=None]
    #right = max(subjct_range)
    #left = min(subjct_range)
    right = max(align_dict)
    left = min(align_dict)

    iLeft = min(align_dict.values())
    iRight = max(align_dict.values())

    if iRight - iLeft + 1 != len(align_dict):
        print('detect deletion!')
        cmn.append_file(name + '\n', 'hasDeletion')

    iLeft = align_dict[left]
    seq = list(seq)
    #print aligned
    while(iLeft != 0):
        iLeft -= 1
        left -= 1
        #print iLeft
        seq[iLeft] = seq[iLeft].lower()
        align_dict[left] = iLeft

    iRight = align_dict[right]
    while(iRight < len(seq) - 1):
        iRight += 1
        right += 1
        seq[iRight] = seq[iRight].lower()
        align_dict[right] = iRight

    aln = []
    for j in range(right):
        #j += 1
        try:
            i = align_dict[j]
            if i == None:
                char = 'N'
            else:
                char = seq[i]
        except KeyError:
            char = '-'
        aln.append(char)
    return ''.join(aln)


def has_indel(pairs):
    #now only discard those of single indel
    for i in range(len(pairs)):
        if i == 0 or i == len(pairs) - 1:
            continue
        aset = [each[0] for each in pairs[i-1:i+2]]
        bset = [each[1] for each in pairs[i-1:i+2]]

        for a0, a1, a2 in [aset, bset]:
            if a1 == None and a0 != None and a2 != None:
                return True
    return False #no indel


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py samlist", file=sys.stderr)
        sys.exit()

    fns = cmn.getid(fn)
    cmn.run('rm hasDeletion 2> /dev/null')

    read_dict = {}
    bad_alignments = []
    seqDict = {}

    for fn in fns:
        print('parsing %s...' % fn)
        try:
            samfile = pysam.AlignmentFile(fn)
        except:
            print('skip empty file %s' % fn)
            continue

        for record in samfile:
            #if record.is_unmapped or record.is_secondary:
            if record.is_unmapped or record.is_secondary:
                continue
            if record.reference_name == '(null)':
                continue

            name = record2name(record)

            #do not allow indel
            #if has_indel(record.aligned_pairs):
            #    print 'skip %s in %s due to indel' % (name, record.reference_name)
            #    continue

            seq = str(record.seq)
            try:
                current = seqDict[name]
                if len(seq) > len(current):
                    seqDict[name] = seq
            except KeyError:
                seqDict[name] = seq


            try:
                read_dict[name].append(record)
            except KeyError:
                read_dict[name] = [record]

    aln_dict = {}
    print('total reads: %s' % len(read_dict))
    for name in read_dict:
        records = read_dict[name]
        seq = records[0].seq
        aligns = [record.aligned_pairs for record in records]
        aligned = consensus_alignment(aligns)
        #print name, aligned
        if aligned == None:
            bad_alignments += records
            print('aligning different regions for %s' % name)
            print('\n'.join(map(str, aligns)))
            continue

        name = '%s(%s)' % (':'.join(name.split(':')[-3:]), len(records))

        aln = seqblock2alignment(aligned, name, seq)
        #print 'align', name, aligned, seq, aln

        aln_dict[name] = aln

    #aln_dict = collapse_same_alignment(aln_dict)
    #output
    #for name in aln_dict:
    #    print 'checkaln_dict', name, aln_dict[name]
    dn = '1st_sam_aln.txt'
    alnDict2output(aln_dict, dn)

    dn = '1st_inconsistent_alignment.txt'
    badDict = {record2name(record, True): seqblock2alignment(record.aligned_pairs, record.query_name, record.seq)
            for record in bad_alignments}
    #print badDict
    alnDict2output(badDict, dn, order='grouping')

    fastas = ['>%s\n%s\n' % (name, seqDict[name]) for name in seqDict]
    dn = '1st_for_denovo.fa'
    cmn.write_file(''.join(fastas), dn)
