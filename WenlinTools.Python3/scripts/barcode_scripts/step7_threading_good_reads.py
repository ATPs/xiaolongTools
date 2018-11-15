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
from collections import Counter
from fullname_lib import get_names_4barcode
gapChars = set(list('-N'))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def merge_into_clusters(clusters, seq):
    seqlen = len(seq)
    keys = list(clusters.keys())
    isMerge = False
    for clID in clusters:
        seq0 = clusters[clID]
        #check difference
        Ndiff = sum([seq0[j] != seq[j] for j in range(seqlen)
            if seq0[j] != '-' and seq[j] != '-'])
        if Ndiff == 0:
            new = []
            isMerge = True
            for j in range(seqlen):
                #a is what is in the current cluster
                a  = seq0[j]
                b = seq[j]
                if a == b:
                    new.append(a)
                elif b == '-':
                    new.append(a)
                elif a == '-':
                    new.append(b)
                else:#different characters
                    print('somthing is wrong!')

            clusters[clID] = ''.join(new)
            #return here to stop the loop
            return isMerge, clusters
    return isMerge, clusters
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_aln(fn):
    seqs = {}
    seq = ''
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue
        name, seq = line.strip().split()
        seqs[name] = seq

    return seqs, len(seq)


def prefix_gap_number(seq):
    i = 0
    while(seq[i] == '-'):
        i += 1
    return i

def aln_hamming_dist(seq1, seq2):
    Ndiff = 0
    overlapN = 0
    for char1, char2 in zip(seq1, seq2):
        if char1 in gapChars or char2 in gapChars:
            continue
        if char1 != char2:
            Ndiff += 1
        overlapN += 1
    return Ndiff, overlapN


def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


def read_and_parse_topHit(fn, refname):
    nameDict = {}
    bad_Ids = []
    good_Ids = []
    bad_baits = {}
    good_baits = {}

    fallbaits = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes.fasta'
    baitSeqs = read_fa(fallbaits)
    for line in cmn.file2lines(fn):
        items = line.split()
        sp, Id, baitID = items[:3]
        print(sp, Id, baitID)
        if refname not in sp:
            bad_Ids.append(Id)
            bad_baits[sp] = baitSeqs[baitID]
        else:
            #refname in name
            if sp.count('|') == 0:
                good_Ids.append(Id)
                good_baits[Id] = baitSeqs[baitID]
        nameDict[Id] = sp
    return nameDict, bad_Ids, good_Ids, good_baits, bad_baits

def bad_Id_number(alist, bad_Ids):
    count = 0
    for name, seq in alist:
        if name in bad_Ids:
            count += 1
    return count

def get_consensus_read(seqs, minCut):
    new = []
    for i in range(len(seqs[0])):
        chars = [seq[i] for seq in seqs
                if seq[i] not in gapChars]
        if len(chars) == 0:
            new.append('-')
            continue
        count_dict = Counter(chars)
        maxChar = max(list(count_dict.keys()), key=lambda x: count_dict[x])
        if count_dict[maxChar] >= minCut:
            new.append(maxChar)
        else:
            new.append('N')
    new = ''.join(new)
    if new.strip('-N') == '':
        return None
    else:
        return new


def prefix_index(seq):
    for i, char in enumerate(seq):
        if char == '-':
            continue
        if char.isupper():
            break
    return i

def suffix_index(seq):
    isStart = False
    for i, char in enumerate(seq):
        if isStart:
            if char == '-' or char.islower():
                break

        if char == '-':
            continue
        if char.isupper():
            isStart = True
    return i

def mapped_percent(seq, primer):
    leftP, rightP = primer
    seq = seq[leftP:658+rightP].strip('-')
    if len(seq) == 0:#not mapped to the barcode region, don't filter out
        return 1
    chars = [char for char in seq
            if char.isupper()]
    fraction = float(len(chars)) / len(seq)

    print(fraction, seq)
    return fraction


def get_concensus(pre, chars, step):
    #print pre, chars
    empty_pre_string = '-' * (step - 1)
    if pre == empty_pre_string:
        pre = None

    #TODO: remove small letter sequence
    if pre != None:
        #filtering reads has a good pre
        pre_chars = []
        for char in chars:
            #currently, if small letter show conflict, didn't consider
            if not isConflict(char[:step-1], pre):
                pre_chars.append(char)
        chars = pre_chars

    good_chars = []
    for char in chars:
        #replace lower character as gap
        newChar = []
        for a in char:
            if a.islower():
                newChar.append('-')
                #newChar.append(a.upper())
            else:
                newChar.append(a)
        good_chars.append(''.join(newChar))

    count_dict = Counter(good_chars)
    empty_string = '-' * step
    try:
        del count_dict[empty_string]
    except:
        pass
    if len(count_dict) == 0:
        return empty_string

    maxChar = max(list(count_dict.keys()), key=lambda x: count_dict[x])
    return maxChar



def number4sorting(seq):
    if not any([char.isupper() for char in seq]):
        return (len(seq), len(seq), len(seq), '')
    i = 0
    i0 = 0
    print(seq)
    while (i < len(seq)):
        if not seq[i].isupper():
            i += 1
            if seq[i] == '-':
                i0 += 1
        else:#is upper now
            break
    j = i
    while (j < len(seq)):
        if seq[i].isupper():
            j += 1
        else:#has lower now
            break
    return (i, j, i0, seq[i:j])

def add_indel(seq, indel_dict):
    seq = list(seq)
    for i, j in indel_dict:
        indel = '-' * len(indel_dict[(i,j)])
        seq[i] += indel
    return ''.join(seq)

def read_indel_info(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        key, a, b, char = line.strip().split('\t')
        i, j = list(map(int, key[1:-1].split(', ')))
        adict[(i,j)] = char
    return adict


def parse_bait(fn):
    alist = []
    adict = {}
    if cmn.filexist('bait_insertion'):
        indel_dict = read_indel_info('bait_insertion')
    else:
        indel_dict = {}

    for line in cmn.file2lines(fn):
        sp, name, seq = line.strip().split()
        if len(indel_dict) != []:
            seq = add_indel(seq, indel_dict)

        adict[name] = seq
        alist.append(name)
    return adict, alist

def search_in_baits(seq, bait_dict, maxN=0):
    #fake it such that only consider mismatch
    seq = seq.replace('-', '').upper()
    names = [name for name in bait_dict
            if seq in bait_dict[name]]
    return names, 0

def cache_search_in_baits(seq, bait_dict, maxN=0):
    global gapChars
    seq = seq.replace('-', '').upper()
    length = len(seq)

    match_dict = {}
    for name in bait_dict:
        baitSeq = bait_dict[name]
        tmp = []
        for i in range(len(baitSeq) - length + 1):
            segment = baitSeq[i:i+length]
            N = sum([segment[j] == seq[j] for j in range(length)
                    if seq[j] not in gapChars])
            tmp.append(N)
        match_dict[name] = min(tmp)

    found = []
    cutN = -1
    #increase allowed mismatches to maxN when nothing was found
    while(cutN <= maxN and len(found) == 0):
        cutN += 1
        names = [name for name in match_dict
                if match_dict[name] <= cutN]
        found += names

    return found, cutN



def isConflict(seq1, seq2):
    length = max(len(seq1), len(seq2))

    for i in range(length):
        try:
            char1 = seq1[i]
            if char1 in gapChars:
                continue
            if char1.islower():
                #treat lower case letter as a gap
                continue

        except IndexError:
            continue

        try:
            char2 = seq2[i]
            if char2 in gapChars:
                continue
            if char2.islower():
                #treat lower case letter as a gap
                continue
        except IndexError:
            continue

        if char1 != char2:
            return True#is bad
    return False

def thread_consensus_seq(good_lines, step=6):
    #step is the length of string considered
    #every time, only move 1 character
    if len(good_lines) == 0:
        return ''
    pre = None
    new = []
    length = max(list(map(len, good_lines)))
    ppSeqs = list(good_lines) # used to removed mark bad sequence
    for i in range(length - step):
        chars = [seq[i:i+step] for seq in ppSeqs
                if len(seq) >= i+step]
        consensus = get_concensus(pre, chars, step)
        #tmp = list([seq for seq in ppSeqs
        #        if not isConflict(seq[i:i+step], consensus)])
        #print 'string:', ''.join(new)
        #print 'cc:', consensus
        #for seq in ppSeqs:
        #    flag = isConflict(seq[i:i+step], consensus)
        #    print 'checking:', i, seq[i:i+step], consensus, flag, seq
        #ppSeqs = list(tmp)
        #print 'ppSeqs', consensus, i, ppSeqs
        pre = consensus[1-step:]
        if new == []:
            new.append(consensus)
        else:
            new.append(consensus[-1])
    return ''.join(new)

def format_name(strF, name):
    if len(name) > 50:
        name = '%s...' % name[:47]

    pname = strF.format(name)
    return pname


def parse_br_name(name_dict, name):
    refname = None
    try:
        defline = name_dict[name]
        if len(defline) > 50:
            if refname in defline:
                defline = '%s|manyOther|%s' % (refname, name)
            else:
                defline = '%s...' % defline[:-3]
        else:
            defline = '%s|%s' % (defline, name)
    except:
        defline = name
    return defline

def get_start_and_end(seq):
    i = 0
    while(not seq[i].isupper()):
        i += 1
    #this i is an upper leter
    j = i
    while(j < len(seq) and seq[j].isupper()):
        j += 1

    return i, j

def isSameSeq(a, b):
    isGood = True
    for char1, char2 in zip(a,b):
        if char1 in gapChars or char2 in gapChars:
            continue
        if char1 != char2:
            isGood = False
            break
    return isGood

def parse_bad_names(name, seq):
    #numbers = list(number4sorting(seq))
    i, j, i0, seq = list(number4sorting(seq))
    items = name.split('|')
    if name[:3] == 'ADD':
        items = items[1:]
    if name[:4] == 'CCrm':
        items = items[1:]
    items = [each for each in items
            if ':' not in each]

    final = [seq, '|'.join(items), i, j, i0]
    #final += numbers
    return final

def spBased_badnames(name, seq):
    sp = name.split('|')[-1].split('(')[0]
    numbers = number4sorting(seq)
    order = list(numbers)
    order.append(sp)
    return tuple(order)


def stack_consensus_seq(seqs):
    if len(seqs) == 0:
        return ''
    length = max(list(map(len, seqs)))
    new = []
    for i in range(length):
        chars = [seq[i] for seq in seqs
                if len(seq) > i and (seq[i] not in gapChars)]
        count_dict = Counter(chars)

        keys = [key for key in count_dict
                if key.isupper()]
        if len(keys) == 0:
            new.append('-')
        else:
            maxChar = max(keys, key=lambda x: count_dict[x])
            sumLength = sum([count_dict[key] for key in keys])
            print('stackRR', i, maxChar, keys, count_dict[maxChar], sumLength)
            if count_dict[maxChar] >= 0.8 * sumLength:
                new.append(maxChar)
            else:
                new.append('N')

    return ''.join(new)


def collapse_same_reads(names, seqDict, hasSp=False):
    spDict = {}

    adict = {}
    for name in names:
        seq = seqDict[name]
        if hasSp:
            sp = name.strip().split('|')[-1].split('(')[0]
            tag = '|' + sp
            name = name.replace(tag, '')
            spDict[name] = sp
        try:
            adict[seq].append(name)
        except KeyError:
            adict[seq] = [name]

    final_dict = {}
    print(spDict)
    for seq in adict:
        alist = adict[seq]
        rpName = alist[0]
        allName = '|'.join(alist)
        if hasSp:
            sp = spDict[rpName]
            rpName = rpName.replace('(', '|%s(' % sp)
            allName = '%s|%s' % (sp, allName)
        print(rpName)
        final_dict[rpName] = allName
    return final_dict


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fbait, fbad, sampleID = sys.argv[1:]
    except:
        print("Usage: *.py filtered_sam_aln.txt fbait bad_alignment.txt", file=sys.stderr)
        sys.exit()

    seqDict, length = read_aln(fn)

    nameDict = get_names_4barcode()
    species = nameDict[sampleID].replace('"', '').replace('\t', '_').replace(' ', '_')

    good_reads = list(seqDict.keys())
    bad_dict, tmp = read_aln(fbad)
    bad_reads = []

    for name in seqDict:
        seq = seqDict[name]
        if not any([char.isupper() for char in seq]):
            bad_dict[name] = seq
            good_reads.remove(name)

    #step4, parse the good reads to be sorted
    good_names = sorted(good_reads, key=lambda x: number4sorting(seqDict[x]))

    #step5, filtering out poorly mapped reads, except it is confident to the bait
    #TODO: currently, those reads are thrown away
    cut_off = 0.25 # requre more than 25% mapped
    good_names = [name for name in good_names
            if (
                #(name in good_Ids) or
                (mapped_percent(seqDict[name], [20,20]) > cut_off)
                )
            ]
    #step5: get the consensus sequence
    good_lines = [seqDict[name] for name in good_names]
    cons_seq = thread_consensus_seq(good_lines, step=6)
    stack_seq = stack_consensus_seq(good_lines)

    #step6: further classify the good reads by the consensus
    bad_cc_names = []
    for name in good_names:
        seq = seqDict[name]
        i, j = get_start_and_end(seq)
        cc_segment = cons_seq[i:j]
        seq = seq[i:j]
        if not isSameSeq(cc_segment, seq):
            bad_cc_names.append(name)

    ### update info
    #for name in bad_cc_names:
    #    good_names.remove(name)
        #try:
        #    pname = name_dict[name]
        #except:
        #    pname = name
        #name_dict[name] = 'CCrm|%s' % name

    #step7, make the output
    final = []
    info = []#report for bad reads
    ###step0, setup the format
    strformat = '{:<54}'

    ###step5.1 get the bait info
    bait_dict, ordered_baits = parse_bait(fbait)
    #bait_dict = read_fa(fbait)
    for name in ordered_baits:
        pname = format_name(strformat, 'b|' + name)
        baitline = '%s%s' % (pname, bait_dict[name])
        final.append(baitline)
        info.append(baitline)

    ###step5.2 format consensus
    name = format_name(strformat, 'threaded_%s' % species)
    ccline = '%s%s' % (name, cons_seq)
    final.append(ccline)
    info.append(ccline)

    name = format_name(strformat, 'stack_%s' % species)
    ccline = '%s%s' % (name, stack_seq)
    final.append(ccline)
    info.append(ccline)

    ###step5.3 (main step): make the lineup
    collapsed_names = collapse_same_reads(good_names, seqDict)
    for name in good_names:
        #Pname = parse_br_name(name_dict, name)
        try:
            Pname = collapsed_names[name]
        except KeyError:
            continue
        Pname = format_name(strformat, Pname)
        line = '%s%s' % (Pname, seqDict[name])
        final.append(line)


    ### report those filtered by bwa mapping to other species
    final.append('#' * 700)
    names = sorted(list(bad_dict.keys()), key=lambda x: spBased_badnames(x, bad_dict[x]))
    collapsed_names = collapse_same_reads(names, bad_dict, True)
    for name in names:
        try:
            Pname = collapsed_names[name]
        except KeyError:
            continue
        Pname = format_name(strformat, Pname)
        line = '%s%s' % (Pname, bad_dict[name])
        final.append(line)

    dn = 'good_read_assembled.txt'
    cmn.write_lines(final, dn)


    #report thoese reads inconsistent with consensus
    for name in bad_cc_names:
        #Pname = parse_br_name(name_dict, name)
        Pname = name
        Pname = format_name(strformat, Pname)
        line = '%s%s' % (Pname, seqDict[name])
        info.append(line)

    info.append('#' * 200)

    #report those reads mapped to other species
    names = sorted(bad_reads, key=lambda x: parse_bad_names(x, seqDict[x]))
    for name in names:
        #Pname = parse_br_name(name_dict, name)
        Pname = name
        Pname = format_name(strformat, Pname)
        line = '%s%s' % (Pname, seqDict[name])
        info.append(line)

    cmn.write_lines(info, 'bad_reads.report')
