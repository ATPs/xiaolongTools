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
        fn, sampleN = sys.argv[1:3]
        sampleN = int(sampleN)
    except:
        print("Usage: *.py fasta sample_size", file=sys.stderr)
        sys.exit()

    seqDict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].split('_')[0]
            else:
                seq = line.strip()
                try:
                    seqDict[name].append(seq)
                except KeyError:
                    seqDict[name] = [seq]

    seqlength = len(seq)

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
    gapcut = 0.8
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

        if N == 0 or N == 1:
            #skip the all gapped positions
            #also skip same character lines
            continue
        else:
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

    indexes = list(range(len(loci_list)))
    goodIndexes = random.sample(indexes, sampleN)

    loci_list = [loci_list[i] for i in goodIndexes]

    new = ['parameters: gapcut %s, Ncut %s' % (gapcut, Ncut)]
    new.append('\n'.join(['loc%s' % each for each in loci_list]))
    new.append('Pop')
    for name in individuals:
        codes = individuals[name]
        codes = [codes[i] for i in goodIndexes]
        line = '%s , %s' % (name, ' '.join(codes))
        new.append(line)

    dn = cmn.lastName(fn).replace('.fa', '').replace('.', '_') + '_genepop_r%s.txt' % sampleN
    new.append('')
    cmn.write_lines(new, dn)

    print('number of loci: %s' % (len(codes)))



