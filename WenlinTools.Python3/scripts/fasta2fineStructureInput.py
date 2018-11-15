#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
import os
current_path = os.path.realpath(__file__)
python_lib = os.path.dirname(current_path)+'/../python_lib'
#print(python_lib)
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_popDef(fpopdef):
    adict = {}
    for fn in cmn.file2lines(fpopdef):
        print(fn)
        popname = cmn.lastName(fn).replace('IDs','')
        IDs = [line.split()[0] for line in cmn.file2lines(fn)]
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

def sortKey(ID):
    global popDict
    try:
        pop = popDict[ID]
    except KeyError:
        pop = -1
    return pop


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fpopDef = sys.argv[1:3]
    except:
        print("Usage: *.py fasta popDef", file=sys.stderr)
        print("popDef is a file containing the files of IDlist", file=sys.stderr)
        sys.exit()

    #according to one paper, 1cM is about 164.6K in butterfly
    #so, 0.01 cM is about 1.6K

    ##############################
    gapCut = 0.1
    Napp = 4
    gapChars = set(list('N-X'))
    ##############################
    outlabel = '%s_gap%s_info%s' % (cmn.lastName(fpopDef), gapCut, Napp)

    #parsing sequence
    seqDict, seqLength = read_fa(fn)
    seqGroups = group_seqDict(seqDict)

    popDict = parse_popDef(fpopDef)

    ID2pop = {}
    for pop in popDict:
        IDs = popDict[pop]
        for ID in IDs:
            ID2pop[ID] = pop

    allkeys = list(seqGroups.keys())
    allkeys = sorted(allkeys, key=lambda x: sortKey(x))

    #parse out id file, also record which is taken
    ids = []
    ordered_keys = []
    Nhap = 0
    for key in allkeys:
        try:
            pop = ID2pop[key]
            label = 1
            Nhap += 2
            ordered_keys.append(key)
        except KeyError:
            pop = 'Ignore'
            label = 0
            continue
        #print key, pop

        ids.append('ind%s %s %s\n' % (key, pop, label))

    dn = outlabel + '.ids'
    cmn.write_file(''.join(ids), dn)


    phase = [str(Nhap)]
    Npos = 0
    Pline = ['P']
    phaseSeqs = [[] for _ in range(Nhap)]

    for i in range(seqLength):
        chars = []
        for key in ordered_keys:
            chars += [seq[i] for seq in seqGroups[key]]

        nonGaps = [char for char in chars
                if char not in gapChars]

        if len(nonGaps) > ((1 - gapCut) * len(chars)):
            Pline.append(i)
            Npos += 1
            for j, char in enumerate(chars):
                if char in gapChars:
                    char = '0'

                phaseSeqs[j].append(char)

    dn = outlabel + '.phase'
    phase.append(str(Npos))
    phase.append(' '.join(map(str, Pline)))
    phase += [''.join(map(str, each)) for each in phaseSeqs]
    phase.append('')
    cmn.write_lines(phase, dn)



