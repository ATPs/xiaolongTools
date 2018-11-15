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
    isGood = True
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
            cutoff = math.ceil(consCut * len(jList))
            if jList.count(checkA) >= cutoff:
                aligned.append((i, checkA))
                alignedN += 1
            else:
                #it is not ok to show conflicts
                aligned.append((i, None))
                conflictN += 1
    print('aligned', alignedN, conflictN, aligned)
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

def seqblock2alignment(aligned, seq):
    #tranform the alinged into a dict for better reference
    align_dict = {}
    for i, j in aligned:
        if j != None:
            align_dict[j] = i

    subjct_range = [i[1] for i in aligned if i[1]!=None]
    right = max(subjct_range)
    left = min(subjct_range)

    iLeft = align_dict[left]
    seq = list(seq)
    #print aligned
    while(iLeft != 0):
        iLeft -= 1
        left -= 1
        #print iLeft
        seq[iLeft] = seq[iLeft].lower()
        align_dict[left] = iLeft


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


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py samlist", file=sys.stderr)
        sys.exit()

    fns = cmn.getid(fn)

    read_dict = {}
    bad_alignments = []
    for fn in fns:
        try:
            samfile = pysam.AlignmentFile(fn)
        except:
            print('skip empty file %s' % fn)
            continue

        for record in samfile:
            if record.is_unmapped or record.is_secondary:
            #if record.is_unmapped:
                continue
            if record.reference_name == '(null)':
                continue

            name = record2name(record, True)

            try:
                read_dict[name].append(record)
            except KeyError:
                read_dict[name] = [record]

    aln_dict = {}
    print('total reads: %s' % len(read_dict))
    for name in read_dict:
        print(name)
        records = read_dict[name]
        seq = records[0].seq
        aligns = [record.aligned_pairs for record in records]
        aligned = consensus_alignment(aligns)
        if aligned == None:
            bad_alignments += records
            print('aligning different regions for %s' % name)
            print('\n'.join(map(str, aligns)))
            continue

        name = '%s(%s)' % (':'.join(name.split(':')[-3:]), len(records))

        aln = seqblock2alignment(aligned, seq)

        aln_dict[name] = aln

    #aln_dict = collapse_same_alignment(aln_dict)
    #output
    dn = '2nd_sam_aln.txt'
    alnDict2output(aln_dict, dn)

    dn = '2nd_inconsistent_alignment.txt'
    badDict = {record2name(record, True): seqblock2alignment(record.aligned_pairs, record.seq)
            for record in bad_alignments}
    print(badDict)
    alnDict2output(badDict, dn, order='grouping')
