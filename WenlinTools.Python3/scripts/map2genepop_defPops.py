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
def chars2line(pos, chars):
    pos = str(pos)
    line = [pos, pos]
    for i in range(0, len(chars), 2):
        a, b = chars[i:i+2]
        if a == b:
            line.append(a)
        else:
            line.append('%s,%s' % (a, b))
    return '\t'.join(line)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def map2seqDict(fn):
    header = []
    seqDict = {}

    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i == 0:
                header = line.strip().split()[2:]
            else:
                chars = line.strip().split()[2:]
                for ii, char in enumerate(chars):
                    name = header[ii]
                    try:
                        seqDict[name].append(char)
                    except KeyError:
                        seqDict[name] = [char]
    seqDict = {name: ''.join(seqDict[name]).replace('N', '-')
            for name in seqDict}
    return seqDict, len(seqDict[name])


def make_code(chars):
    global missing_data, char_label, current_count

    if chars[0] in missing_data:
        code = '0000'
    else:#not gap
        code = []
        for char in chars:
            try:
                subCode = char_label[char]
            except KeyError:
                subCode = str(current_count).zfill(2)
                char_label[char] = subCode
                current_count += 1
            code.append(subCode)
        code = ''.join(code)
    return code


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py map pop1 pop2 ...", file=sys.stderr)
        sys.exit()

    poplist = sys.argv[2:]

    if len(poplist) == 0:
        print('please specify populations!')
        sys.exit()

    popdict = {}
    for fpop in poplist:
        IDs = cmn.getid(fpop)
        name = cmn.lastName(fpop)
        popdict[name] = IDs
        #for ID in IDs:
        #    popdict[ID] = name

    #seqDict = {}
    #with open(fn) as fp:
    #    for line in fp:
    #        if line[0] == '>':
    #            name = line[1:].split('_')[0]
    #        else:
    #            seq = line.strip()
    #            try:
    #                seqDict[name].append(seq)
    #            except KeyError:
    #                seqDict[name] = [seq]

    #seqlength = len(seq)
    seqDict, seqlength = map2seqDict(fn)

    names = list(seqDict.keys())
    char_label = {}
    current_count = 1
    missing_data = set('N -'.split())
    #header = ['fake', 'position'] + names
    #header = '\t'.join(header)

    seqs = []
    #hetero = [[header], [header], [header], [header]]
    for name in names:
        seqs += seqDict[name]

    loci_list = []
    individuals = {}
    gapcut = 0.3
    #enrich more informative positions
    #each nucleatide appear should have larger count than this
    #as this code was designed for one population, so just ignore this
    Ncut = 1
    for i in range(len(seqs[0])):
        chars = [seq[i] for seq in seqs]
        check_chars = [char for char in chars if (char != '-') and (char != 'N')]
        if len(check_chars) < (gapcut * len(chars)):
            continue
        N = len(set(check_chars))
        count_dict = Counter(check_chars)
        if any([countN < Ncut for countN in list(count_dict.values())]):
            continue

        #if N == 0 or N == 1:
            #skip the all gapped positions
            #also skip same character lines
        #    continue
        #else:
        if True:
            #output
            loci_list.append(i)
            for name in seqDict:
                seqs = seqDict[name]
                chars = [seq[i] for seq in seqs]
                code = make_code(chars)
                try:
                    individuals[name].append(code)
                except KeyError:
                    individuals[name] = [code]

    print(list(individuals.keys()))

    new = ['parameters: gapcut %s, Ncut %s' % (gapcut, Ncut)]
    new.append('\n'.join(['loc%s' % each for each in loci_list]))
    for pop in popdict:
        samples = popdict[pop]
        new.append('Pop')
        for name in samples:
            Sname = '%s_%s' % (pop, name)
            line = '%s , %s' % (Sname, ' '.join(individuals[name]))
            new.append(line)

    dn = cmn.lastName(fn).replace('.fa', '').replace('.', '_') + '_genepop.txt'
    new.append('')
    cmn.write_lines(new, dn)

    print('number of loci: %s' % (len(individuals[name]) ))



