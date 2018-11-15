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
gapChars = set('X - N'.split())

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def merge_into_clusters(clusters, seq):
    seqlen = len(seq)
    #keys = clusters.keys()
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


def read_and_parse_topHit(fn, refname):
    nameDict = {}
    bad_Ids = []
    good_Ids = []
    bad_baits = {}
    good_baits = {}

    fallbaits = '/work/00412/mtang/sequencing/scripts/data/barcodes/all_barcodes.fasta'
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


def number4sorting(seq):
    if seq.islower():
        return (len(seq), len(seq), len(seq), '')
    i = 0
    i0 = 0
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

def parse_bait(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        sp, name, seq = line.strip().split()
        adict[name] = seq
    return adict

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


def format_name(strF, name, indel_list=[]):
    global hasDelLabel
    if len(indel_list) != 0:
        labels = set(name.split('|'))
        if len(labels & indel_list) != 0:
            name = '[hasDel]' + name
            hasDelLabel = True
    if len(name) > 50:
        name = '%s...' % name[:47]

    pname = strF.format(name)
    return pname

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict

def read_baits(fn):
    adict = {}
    hasPrimer = True
    new = []
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue
        sp, name, seq = line.split()
        if len(seq) != 698:
            hasPrimer = False
            if len(seq) == 658:
                #fixable
                seq = add_primer(seq)
            else:
                print('Error! didn\'t recognize the length of the bait %s %s' % (sp, name))
                sys.exist()
        newline = '%s\t%s\t%s\n' % (sp, name, seq)
        new.append(newline)
        key = '%s_%s' % (sp, name)
        adict[key] = seq

    if not hasPrimer:
        print('revise the input baits to add primer...')
        cmn.write_file(''.join(new), fn)

    return adict


def add_primer(seq):
    return 'ACTAATCATAAAGATATTGG%sTGATTTTTTGGTCATCCAGA' % seq.strip()

def parse_ref(seqDict):
    cmn.mkdir('baits')

    newDict = {}
    for i, name in enumerate(seqDict):
        seq = seqDict[name]
        fnlabel = 'bait%s' % i
        dn = 'baits/%s.fa' % fnlabel
        name = name.replace('*', '').replace('"', "'")
        fasta = '>%s\n%s\n' % (name, seq)
        cmn.write_file(fasta, dn)
        cmd = 'module add bwa; bwa index %s -p %s' % (dn, fnlabel)
        cmn.run(cmd)
        newDict[name] = dn
    return newDict

def group_fq(fqlist):
    adict = {}
    for fn in fqlist:
        fnlabel = cmn.lastName(fn)
        sp = fnlabel.split('_')[0]
        if sp not in adict:
            adict[sp] = [None, None, None]
        if '_R1' in fnlabel:
            adict[sp][0] = fn
        elif '_R2' in fnlabel:
            adict[sp][1] = fn
        elif '_singleton' in fnlabel:
            adict[sp][2] = fn
        else:
            print('Error! Can not recognize file %s' % fn)
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
            #print 'stackRR', i, maxChar, keys, count_dict[maxChar], sumLength
            if count_dict[maxChar] >= 0.8 * sumLength:
                new.append(maxChar)
            else:
                new.append('N')

    return ''.join(new)


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


def read_assembled_file(fn):
    baits = {}
    ordered_baits = []
    good_reads = {}
    junk_reads = {}
    isJunk = False
    sampleLabel = ''

    for line in cmn.file2lines(fn):
        items = line.split()
        try:
            name, seq = items
        except:
            if '######' in line:
                isJunk = True
                continue

        if name[:2] == 'b|':
            baits[name] = seq
            ordered_baits.append(name)
        elif name.startswith('stack_'):
            stack_seq = seq

        elif name.startswith('threaded_'):
            thread_seq = seq

            try:
                #sampleID = line.strip().split()[0].split('_')[2]
                sampleLabel = '_'.join(line.strip().split()[0].split('_')[1:])
            except:
                pass

        else:
            if isJunk:
                junk_reads[name] = seq
            else:
                #not reaching the junk bin yet
                #filter out those with no more than 20 letters
                Nup = len([each for each in seq if each.isupper()])
                if Nup > 20:
                    good_reads[name] = seq
                else:
                    junk_reads[name] = seq

    return baits, ordered_baits, stack_seq, thread_seq, good_reads, junk_reads, sampleLabel


def rescue_junk_reads(ccSeq, junk_reads):
    still_junk = {}
    rescued = {}
    for name in junk_reads:
        seq = junk_reads[name]
        Nup = len([each for each in seq if each.isupper()])
        if isConflict(seq, ccSeq) or Nup <= 20:
            still_junk[name] = seq
        else:
            rescued[name] = seq
    return rescued, still_junk

def isConflict(seq1, seq2):
    #seq2 is the ccSeq
    global misN

    length = max(len(seq1), len(seq2))

    wrongN = 0
    gapN = 0
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
                gapN += 1
                continue
            if char2.islower():
                gapN += 1
                #treat lower case letter as a gap
                continue
        except IndexError:
            continue

        if char1 != char2:
            wrongN += 1

    length = len(seq1.strip('-'))
    print('checkLog', gapN, length, seq1)
    if gapN > length * 0.2:
        print('didn\'t rescue the read covering the gap region: %s' % seq1)
        return True
    else:
        if wrongN > misN:#make it such when misN=0 means no mismatch allowed
            return True#is bad
        else:
            return False


def filter_clean_reads(readDict):
    rdict = {}
    for name in readDict:
        seq = readDict[name]
        checkSeq = seq.strip('N- ')
        N = len(checkSeq)
        mapN = sum([each.isupper() for each in checkSeq])
        if N >= 40 and float(mapN) / N >= 0.5:
            rdict[name] = seq
    return rdict


def isGoodEnd(seq, side):
    if side == 'left':
        checkSeq = seq[:20].strip('-')
    elif side == 'right':
        checkSeq = seq[678:].strip('-')
    if len(checkSeq) == 0:
        #no left sequence, then it is good
        return True
    else:
        #has left sequence, need to check aligned fraction
        Ntotal = check_aligned_seq(seq)
        Nleft = check_aligned_seq(checkSeq)
        if Nleft / float(Ntotal) > 0.5:
            return False
        else:
            return True


def check_aligned_seq(seq):
    N = 0
    for char in seq:
        if char in gapChars:
            continue
        if char.isupper():
            N += 1
    return N


def rejecting_badPrimer_sequence(good_reads):
    stillGood, badEnds = {}, {}
    for name in good_reads:
        seq = good_reads[name]
        isGood = True

        isGood = isGood & isGoodEnd(seq, 'left')
        isGood = isGood & isGoodEnd(seq, 'right')

        if isGood:
            stillGood[name] = seq
        else:
            print('rejecting this end read: %s' % name)
            badEnds[name] = seq

    return stillGood, badEnds



if __name__=='__main__':
    #options=parse_options()
    try:
        fn =sys.argv[1]
    except:
        print("Usage: *.py good_reads_assembled [allowed_mismatch=0]", file=sys.stderr)
        sys.exit()

    try:
        misN = int(sys.argv[2])
    except:
        misN = 1


    hasDelLabel = False
    #NOTE: new feature: rejecting reads that matched mostly to the extended ends
    #NOTE: change it such that bad reads are not rescued due to gaps
    if cmn.filexist('hasDeletion'):
        indel_list = set(cmn.file2lines('hasDeletion'))
    else:
        indel_list = set([])

    #add primer if not added
    bait_dict, ordered_baits, stack_seq, thread_seq, good_reads, junk_reads, sampleLabel = read_assembled_file(fn)
    print('junk1', len(junk_reads))

    #use stack to rescue
    all_rescued = {}
    while(True):
        #this means more reads are rescued
        #continue doing rescue

        rescued_reads, junk_reads = rescue_junk_reads(stack_seq, junk_reads)
        if len(rescued_reads) == 0:
            break

        print('loop junk', len(junk_reads))
        all_rescued.update(rescued_reads)
        #recompute stack consensus
        stack_seq = stack_consensus_seq(list(all_rescued.values()) + list(good_reads.values()))

    good_reads.update(all_rescued)
    good_reads, end_reads = rejecting_badPrimer_sequence(good_reads)
    junk_reads.update(end_reads)

    stack_seq = stack_consensus_seq(list(good_reads.values()))
    thread_seq = thread_consensus_seq(list(good_reads.values()))

    clean_reads = filter_clean_reads(good_reads)
    clean_seq = stack_consensus_seq(list(clean_reads.values()))

    bad_dict = junk_reads
    ###step0, setup the format
    final = []
    info = []
    strformat = '{:<54}'

    ###step5.1 get the bait info
    for name in ordered_baits:
        pname = name
        #if len(indel_list) != 0:
        if cmn.filexist('bait_insertion'):
            pname += '[hasInsert]'
        pname = format_name(strformat, pname)
        baitline = '%s%s' % (pname, bait_dict[name])
        final.append(baitline)
        info.append(baitline)

    #sampleLabel = ''
    #if sampleID != '':
    #    sampleLabel = '_' + sampleID

    ###step5.2 format consensus
    name = format_name(strformat, 'threaded_%s' % sampleLabel)
    ccline = '%s%s' % (name, thread_seq)
    final.append(ccline)
    info.append(ccline)

    name = format_name(strformat, 'stack_%s' % sampleLabel)
    ccline = '%s%s' % (name, stack_seq)
    final.append(ccline)
    info.append(ccline)

    name = format_name(strformat, 'clean_%s' % sampleLabel)
    ccline = '%s%s' % (name, clean_seq)
    final.append(ccline)
    info.append(ccline)

    ###step5.3 (main step): make the lineup
    #collapsed_names = collapse_same_reads(good_names, seqDict)
    good_names = sorted(good_reads, key=lambda x: number4sorting(good_reads[x]))
    collapsed_names = good_reads
    for name in good_names:
        #Pname = parse_br_name(name_dict, name)
        #try:
        #    Pname = collapsed_names[name]
        #except KeyError:
        #    continue
        Pname = name
        Pname = format_name(strformat, Pname, indel_list)
        line = '%s%s' % (Pname, good_reads[name])
        final.append(line)


    ### report those filtered by bwa mapping to other species
    final.append('#' * 700)
    names = sorted(list(bad_dict.keys()), key=lambda x: spBased_badnames(x, bad_dict[x]))
    #collapsed_names = collapse_same_reads(names, bad_dict, True)
    for name in names:
        #try:
        #    Pname = collapsed_names[name]
        #except KeyError:
        #    continue
        Pname = format_name(strformat, name, indel_list)
        line = '%s%s' % (Pname, bad_dict[name])
        final.append(line)

    if hasDelLabel:
        dn = 'rescued_read_assembled_mis%s_withDeletion.txt' % misN
    else:
        dn = 'rescued_read_assembled_mis%s.txt' % misN

    cmn.write_lines(final, dn)

    Ngood = len(good_reads)
    Njunk = len(junk_reads)
    #print junk_reads
    statInfo = 'junk:good = %s:%s(%s)\n' % (Njunk, Ngood, float(Njunk)/(Njunk + Ngood))
    cmn.write_file(statInfo, 'rescued_ratio_mis%s.txt' % misN)

