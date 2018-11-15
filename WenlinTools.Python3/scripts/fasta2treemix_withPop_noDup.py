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


def parse_popDef(fpopdef):
    adict = {}
    for fn in cmn.file2lines(fpopdef):
        print(fn)
        popname = cmn.lastName(fn).replace('IDs','')
        IDs = cmn.file2lines(fn)
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fpopdef = sys.argv[1:]
    except:
        print("Usage: *.py fasta popdef", file=sys.stderr)
        sys.exit()


    seqDict = read_fa(fn)

    popDict = parse_popDef(fpopdef)

    SNPdict = {}

    for name in seqDict:
        pop = tell_pop(name, popDict)
        if pop == None:
            print('%s is not included in population definition, skip' % name)
            continue

        seq = seqDict[name]

        if pop not in SNPdict:
            SNPdict[pop] = []

        for i, char in enumerate(seq):
            try:
                SNPdict[pop][i].append(char)
            except IndexError:
                SNPdict[pop].append([char])

    length = len(seq)

    header = list(SNPdict.keys())

    new = ['\t'.join(header)]

    for i in range(length):
        popSNPs = [SNPdict[pop][i] for pop in header]
        allSNPs = sum(popSNPs, [])
        count_dict = Counter(allSNPs)
        for gapChar in list('N-'):
            try:
                del count_dict[gapChar]
            except:
                pass

        #decide if it is bialletic
        if len(count_dict) == 2:
            keys = list(count_dict.keys())
            line = []
            for popSNP in popSNPs:
                counts = [ str(popSNP.count(key))
                        for key in keys]
                line.append(','.join(counts))

            new.append('\t'.join(line))

    new.append('')
    dn = fn.replace('.fa', '') + '_treemix.snp'
    cmn.write_lines(new, dn)



