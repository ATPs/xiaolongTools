#!/usr/bin/env python

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
from collections import Counter
gapChars = set(['X', '-', 'N'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_letter_file(fn):
    #3614_mito       63      A       T       (phase:63)
    #3614_mito       64      A       T       (phase:63, PhasingInconsistent)
    adict = {}
    bad_p = set([])
    for line in cmn.file2lines(fn):
        items = line.strip().split('\t')
        scaf, index, char1, char2, phaseInfo = items
        index = int(index)

        if scaf not in adict:
            adict[scaf] = {}

        if 'PhasingInconsistent' in phaseInfo:
            bad_p.add((scaf, index))
            phaseInfo = phaseInfo.replace(', PhasingInconsistent', '')

        phase = phaseInfo[1:-1].split(':')[-1]
        adict[scaf][index] = (char1, char2, phase)
    return adict, bad_p

def determine_contigs(adict, bad_p):
    contigs = {}
    ordered_list = []
    for scaf in adict:
        subdict = adict[scaf]
        count = 1
        keys = list(subdict.keys())
        keys.sort()
        current_contigs = []
        for index in keys:
            char1, char2, phase = subdict[index]

            if phase == 'inconsistent':
                #found inconsistent phase, skip this positions
                #but if there is anything in current_contigs, output it
                if len(current_contigs) != 0:
                    contig_key = '%s_%s' % (scaf, count)
                    count += 1
                    ordered_list.append(contig_key)
                    contigs[contig_key] = list(current_contigs)
                    current_contigs = []
                continue


            if len(current_contigs) == 0:
                current_contigs.append((scaf, index, char1, char2, phase))
                continue
            else:
                isBreak = False
                if (scaf, index) in bad_p:
                    isBreak = True
                else:
                    current_phase = current_contigs[0][-1]
                    if phase != current_phase:
                        isBreak = True

                if isBreak:
                    contig_key = '%s_%s' % (scaf, count)
                    count += 1
                    ordered_list.append(contig_key)
                    contigs[contig_key] = list(current_contigs)
                    current_contigs = [(scaf, index, char1, char2, phase)]
                else:
                    current_contigs.append((scaf, index, char1, char2, phase))
        if len(current_contigs) != 0:
            contig_key = '%s_%s' % (scaf, count)
            ordered_list.append(contig_key)
            contigs[contig_key] = list(current_contigs)
    return contigs, ordered_list
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_samfile(fn):
    samfile = pysam.AlignmentFile(fn)

    rdict = {}
    for record in samfile:
        name = record.query_name
        if record.is_unmapped or record.is_secondary:
            continue

        aligned = record.get_aligned_pairs(matches_only=True)

        length = len(record.query_sequence)

        #assuming new sample, so remove short reads
        if length < 50:
            continue

        if len(aligned) < 0.5 * length:
            continue

        try:
            alist = rdict[name]
        except KeyError:
            alist = [None, None]

        if record.is_read1:
            alist[0] = record
        elif record.is_read2:
            alist[1] = record
        else:
            alist[0] = record

        rdict[name] = alist
    return rdict


def tell_consistent_phases(sorted_list):
    #sorted list from smallest to largest
    newlist = []
    for each in sorted_list[::-1]:
        if len(newlist) == 0:
            newlist.append(each)
        else:
            conflict_flag = False
            for phase in newlist:
                if len(set(phase) & set(each)) != 0:
                    conflict_flag = True
                    break

            if not conflict_flag:
                newlist.append(each)
    return newlist


def check_linked_chars(paired_samDict, p1, p2, scaf):
    global cons_seq
    #TODO: use scaf
    #I think positions in p1 are always smaller than p2
    linked_p = {}
    p1 = set(p1)
    p2 = set(p2)
    covered_count = 0
    for name in paired_samDict:
        read1, read2 = paired_samDict[name]
        #TODO: I think they use python index, in order to compare, I need to shift them
        #{aligned_position: char}
        aligned_dict = parse_aligned_reads([read1, read2])
        keys = set(aligned_dict.keys())
        #print keys, p1, p2
        covered_p1 = keys & p1
        covered_p2 = keys & p2
        if len(covered_p1) != 0 and len(covered_p2) != 0:
            covered_count += 1
            for pos1 in covered_p1:
                char1 = aligned_dict[pos1]
                for pos2 in covered_p2:
                    char2 = aligned_dict[pos2]
                    key = (pos1, pos2)
                    value = (char1, char2, name)
                    try:
                        linked_p[key].append(value)
                    except KeyError:
                        linked_p[key] = [value]

    #TODO:filtering the linked positions if they are not bi
    new = {}
    for key in linked_p:
        #firstly, filtering to get only the top 2 abundant sequence
        alist = linked_p[key]
        if len(alist) == 0:
            continue

        blist = [each[:2] for each in alist]
        count_dict = Counter(blist)
        sort_list = sorted(list(count_dict.keys()), key=lambda x: count_dict[x])

        if len(count_dict) == 1:
            new[key] = alist
            continue

        #from the most to least, find the consistent phases
        consistent_phases = tell_consistent_phases(sort_list)

        phases = sort_list[-2:]

        if sum([count_dict[phase] for phase in phases]) > (0.8 * len(blist)):
            newlist = []
            for phase in consistent_phases:
                newlist += [each for each in alist
                        if each[:2] == phase]
            new[key] = newlist
            continue

        #reach here if the most dominiated two copies didn't make 80% of the phase
        #usually, the number of linked type is larger than two
        #then we need to look into the reads and take the ones most similar to query
        diff_dict = {}
        for char1, char2, name in alist:
            reads = paired_samDict[name]
            Ndiff = tell_diff_chars(reads, cons_seq)
            try:
                diff_dict[Ndiff].append(name)
            except KeyError:
                diff_dict[Ndiff] = [name]

        Nmin = min(diff_dict.keys())
        good_names = set(diff_dict[Nmin])
        print('NminReport: filtered %s from %s reads' % (len(good_names), len(alist)))

        newlist = [each for each in alist
                if alist[-1] in good_names]
        new[key] = newlist
        #check the sequence identity of the read to the consensus sequence
        #only take the top identity group
        #if this removes too much, report a warnning
    return new

def tell_diff_chars(reads, cons_seq):

    Ndiff = 0
    for read in reads:
        if read == None:
            continue

        scaf = read.reference_name
        ccSeq = cons_seq[scaf]
        seq = read.query_sequence
        aligned = read.get_aligned_pairs(matches_only=True)
        #print 'comparing read: %s %s' % (read.query_name, ''.join([seq[each[0]] for each in aligned]))
        #print 'comparing ref: %s' % (''.join([ccSeq[each[1]] for each in aligned]))
        tmp = 0
        for i, j in aligned:
            char1 = seq[i]
            char2 = ccSeq[j]
            if char2 in gapChars:
                continue
            else:
                if char1 != char2:
                    Ndiff += 1
                    tmp += 1
        #print 'diff positions: %s' % tmp
    return Ndiff


def parse_aligned_reads(alist):
    adict = {}
    for each in alist:
        if each == None:
            continue
        aligned = each.get_aligned_pairs(matches_only=True)
        seq = each.query_sequence

        for i, j in aligned:
            char = seq[i]
            pos = j + 1 #TODO:needed or not?
            adict[pos] = char
    return adict


def check_support_positions(alist, blist, positions):
    N = 0
    for key in positions:
        char_list = positions[key]
        #for every position pair, the total weight is 1
        suppN = 0
        for chars in char_list:
            isSupport1 = tell_ifSupport(alist, key, chars)
            isSupport2 = tell_ifSupport(blist, key, chars)
            if isSupport1 or isSupport2:
                suppN += 1

        if len(char_list) == 0:
            continue

        N += float(suppN) / len(char_list)
    return N


def tell_ifSupport(alist, key, chars):
    pos1, pos2 = key
    ref_char1, ref_char2, name = chars
    obs_char1 = alist[0][pos1]
    obs_char2 = alist[1][pos2]
    if ref_char1 == obs_char1 and (ref_char2 == obs_char2):
        return True
    else:
        return False


def generate_two_phases(contig):
    ph1, ph2 = {}, {}
    for each in contig:
        scaf, index, char1, char2, phase = each
        ph1[index] = char1
        ph2[index] = char2
    return ph1, ph2

def test_possible_recombination(linked_positions, contig1, contig2):
    #return values includes:
    #1. swap: swap current phase
    #2. keep: keep current phase
    #3. unknown: not sure
    phase11, phase12 = generate_two_phases(contig1)
    phase21, phase22 = generate_two_phases(contig2)

    support1 = check_support_positions([phase11, phase21], [phase12, phase22],  linked_positions)
    support2 = check_support_positions([phase11, phase22], [phase12, phase21],  linked_positions)

    cutoff = 0.8
    try:
        noSwapP = float(support1) / (support1 + support2)
    except ZeroDivisionError:
        return 'unknown'
    #SwapP = float(support1) / len(linked_positions)
    print(len(linked_positions), support1, support2, noSwapP)
    if noSwapP > cutoff:
        return 'keep'
    elif noSwapP < (1 - cutoff):
        return 'swap'
    else:
        return 'unknown'


def correct_false_snp_call(letter_dict, covDict,inconsistent_positions):
    cutoff = 0.9

    new = {}
    corrected = {}
    for scaf in letter_dict:
        new[scaf] = {}
        corrected[scaf] = {}
        pdict = letter_dict[scaf]
        for pos in pdict:
            char1, char2, phase = pdict[pos]
            try:
                cov1 = covDict[scaf][pos][char1]
            except KeyError:
                cov1 = 0
            try:
                cov2 = covDict[scaf][pos][char2]
            except KeyError:
                cov2 = 0
            total = float(cov1 + cov2)

            if (scaf, pos) in inconsistent_positions:
                if cov1 > cov2:
                    char = char1
                    cov = cov1
                else:
                    cov = cov2
                    char = char2
                corrected[scaf][pos] = (char, cov, '(forced) %s' % total)
                new[scaf][pos] = (char1, char2, 'inconsistent')

            elif total == 0:
                #no high quality read covers it, skip this position
                print('WARNNING: no high quality read covered for %s %s, drop it' % (pos, pdict[pos]))
                corrected[scaf][pos] = ('N', cov1, total)

            elif cov1/total > cutoff:
                corrected[scaf][pos] = (char1, cov1, total)
            elif cov2/total > cutoff:
                corrected[scaf][pos] = (char2, cov2, total)
            else:#no need to correct
                new[scaf][pos] = (char1, char2, phase)
    return new, corrected



def compute_coverage_from_sam(paired_samDict):
    covDict = {}
    for name in paired_samDict:
        reads = paired_samDict[name]

        for read in reads:
            if read == None:
                continue

            scaf = read.reference_name
            if scaf not in covDict:
                covDict[scaf] = {}

            seq = read.query_sequence
            aligned = read.get_aligned_pairs(matches_only=True)
            for i, j in aligned:
                index = j + 1
                char = seq[i]
                if index not in covDict[scaf]:
                    covDict[scaf][index] = {}

                try:
                    covDict[scaf][index][char] += 1
                except KeyError:
                    covDict[scaf][index][char] = 1

    return covDict

def make_mito_from_covDict(covDict, maxIndex):
    new = {}
    for scaf in covDict:
        scaf_cov = covDict[scaf]
        #maxIndex = max(scaf_cov.keys())
        seq = []
        for i in range(maxIndex):
            try:
                cov = scaf_cov[i+1]
            except KeyError:
                print('warnning: missing coverage for index %s' % (i+1))
                seq.append('N')
                continue

            maxChar = max(list(cov.keys()), key=lambda x: cov[x])
            #if cov[maxChar] > (0.9 * sum(cov.values())):
            #    seq.append(maxChar)
            #else:
            #    seq.append('N')
            seq.append(maxChar)

        new[scaf] = ''.join(seq)
    return new

def filtered_chars_near_gap(linked_p, paired_samDict):
    new = {}
    for key in linked_p:
        alist = linked_p[key]
        positions = set(key)
        newlist = []
        for char1, char2, name in alist:
            reads = paired_samDict[name]
            isGood = True
            for read in reads:
                if read == None:
                    continue
                aligned = read.aligned_pairs
                aligned_dict = {j+1: i for i, j in aligned
                        if j != None}
                half_window = 2

                for pos in positions:
                    try:
                        i = aligned_dict[pos]
                    except KeyError:
                        continue

                    if i == None:
                        isGood = False
                    else:
                        testSet = []
                        for jj in range(pos-half_window,pos+half_window+1):
                            try:
                                testSet.append(aligned_dict[jj])
                            except KeyError:
                                continue

                        if None in testSet:
                            isGood = False
            if isGood:
                newlist.append((char1, char2, name))
            else:
                print('filtered because of gap: %s %s' % (key, name))
        new[key] = newlist
    return new


def swap_phases(contig):
    #contig_info = [(scaf, position1, A, T, phase), ()]
    newlist = []
    for each in contig:
        scaf, position1, char1, char2, phase = each
        phase += '[swap]'
        newlist.append((scaf, position1, char2, char1, phase))
    return newlist


if __name__=='__main__':
    #options=parse_options()
    try:
        fsam, mitoLength = sys.argv[1:]
        maxN = int(mitoLength)
    except:
        print("Usage: *.py *.sam mitoLength", file=sys.stderr)
        sys.exit()

    outlabel = cmn.lastName(fsam)[:-4]
    print(outlabel)
    #{read_query_name: [record1, record2]}
    paired_samDict = read_samfile(fsam)
    #covDict[scaf][index][char]
    covDict = compute_coverage_from_sam(paired_samDict)
    cons_seq = make_mito_from_covDict(covDict, maxN)
    seq = []
    for scaf in cons_seq:
        seq.append(cons_seq[scaf])

    seq = ''.join(seq)

    dn = '%s_mito.fa' % outlabel
    fasta = '>%s_mito\n%s\n' % (outlabel, seq)
    cmn.write_file(fasta, dn)
