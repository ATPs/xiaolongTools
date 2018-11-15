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
from collections import Counter
import random


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_popDef(fpopdef, freceipt, inclusion):
    adict = {}
    receipt = cmn.lastName(freceipt)
    for fn in cmn.file2lines(fpopdef):
        print(fn)
        if cmn.lastName(fn) == receipt:
            continue
        popname = cmn.lastName(fn).replace('IDs','').rstrip('_')
        IDs = [line.split()[0] for line in cmn.file2lines(fn)]
        IDs = set(IDs) & inclusion
        adict[popname] = IDs
    return adict

def tell_pop(name, popDict):
    pop = None
    ID = name.split('_')[0]
    for each in popDict:
        IDs = set(popDict[each])
        if ID in IDs:
            pop = each
            break

    #if pop == None:
    #    print 'Error! can not find pop for %s' % name
    #    sys.exit()
    #else:
    #    return pop
    return pop


def group_seqDict(adict):
    rdict = {}
    for name in adict:
        key = name.split('_')[0]
        seq = adict[name]
        try:
            rdict[key].append(seq)
        except:
            rdict[key] = [seq]
    return rdict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict, len(seq)


def take_position_chars(seqGroups, keys, i):
    chars = []
    for key in keys:
        chars += [seq[i] for seq in seqGroups[key]]
    return chars

if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fdonor, freceipt = sys.argv[1:4]
    except:
        print("Usage: *.py fasta donorPopDef receiptIDs", file=sys.stderr)
        sys.exit()

    #so 10M should be 1 morgan

    ##############################
    gapCut = 0.1
    Napp = 4
    gapChars = set(list('N-X'))
    linkage = 10000000.0 #the larger the stronger linkage
    ##############################
    outlabel = 'chp_%s_%s_%s_gap%s_info%s' % (cmn.lastName(fn), cmn.lastName(fdonor), cmn.lastName(freceipt), gapCut, Napp)

    #parsing sequence
    seqDict, seqLength = read_fa(fn)
    seqGroups = group_seqDict(seqDict)
    included_IDs = set(seqGroups.keys())

    donor_dict = parse_popDef(fdonor, freceipt, included_IDs)
    donor_keys = list(donor_dict.keys())
    random.shuffle(donor_keys)
    print(donor_keys)
    donorIDs = []

    donorF = []
    key_groups = [] #this is used to fill in gap for each group
    for key in donor_keys:#just to garantee ordering
        IDs = donor_dict[key]
        key_groups.append(IDs)
        donorIDs += list(IDs)
        line = 'p%s %s\n' % (key, len(IDs)*2)
        donorF.append(line)

    dn = outlabel + '.donor'
    cmn.write_file(''.join(donorF), dn)

    receiptIDs = [line.split()[0] for line in cmn.file2lines(freceipt)]
    receiptIDs = list(set(receiptIDs) & included_IDs)
    dn = outlabel + '_ind_record.list'
    cmn.write_lines(receiptIDs, dn)
    key_groups.append(receiptIDs)

    ids = []
    ordered_keys = donorIDs + receiptIDs

    phase = [str(len(donorIDs)*2)]
    phase.append(str(len(donorIDs) + len(receiptIDs)))

    Npos = 0
    Pline = ['P']
    phaseSeqs = [[] for _ in range(len(ordered_keys)*2)]

    for i in range(seqLength):
        chars = []
        isBad = False
        for keys in key_groups:
            subchars = take_position_chars(seqGroups, keys, i)

            #require half of subchars is not gap
            nonGap_sub = [char for char in subchars
                    if char not in gapChars]
            #print i, subchars, nonGap_sub
            if len(nonGap_sub) < len(subchars) * 0.5:
                isBad = True
                break

            #random sample the positions to fill in gap
            newSubChars = []
            for char in subchars:
                if char in gapChars:
                    newChar = random.sample(nonGap_sub, 1)[0]
                else:
                    newChar = char
                newSubChars.append(newChar)

            chars += newSubChars
        #print 'chars', chars
        if isBad:
            continue
        nonGaps = [char for char in chars
                if char not in gapChars]
        #print i, chars
        #print i, nonGaps
        if len(nonGaps) > ((1 - gapCut) * len(chars)):
        #if len(nonGaps) == len(chars):
            #filtering for infoP
            count_dict = Counter(nonGaps)
            if len(count_dict) != 2:
                continue
            if any([count_dict[key] < Napp for key in count_dict]):
                continue
            #print i, 'isGood'
            Pline.append(i)
            Npos += 1
            for j, char in enumerate(chars):
                if char in gapChars:
                    char = '0'

                phaseSeqs[j].append(char)

    dn = outlabel + '.hap'
    phase.append(str(Npos))
    phase.append(' '.join(map(str, Pline)))
    print('number of positions: %s' % (len(Pline) -1))
    phase.append('S' * Npos)
    phase += [''.join(map(str, each)) for each in phaseSeqs]
    phase.append('')
    cmn.write_lines(phase, dn)

    #parse out recomb
    linkage = float(linkage)
    positions = Pline[1:]
    recomb = ['pos morgan.dist']
    for i in range(len(positions) - 1):
        p1 = positions[i]
        p2 = positions[i+1]
        dist = p2 - p1
        morgan = dist / linkage
        recomb.append('%s %s' % (p1, morgan))
    recomb.append('%s 0' % (positions[-1]))
    recomb.append('')
    dn = outlabel + '.recomb'
    cmn.write_lines(recomb, dn)

