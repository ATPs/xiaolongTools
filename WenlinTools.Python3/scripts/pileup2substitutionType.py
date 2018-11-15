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


reverse_chars = {
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'T': 'A'
        }
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    #currently, deal with three scenario
    #1. insertion and deletion: show sugguested change
    #2. heterozyocity position: show hetero ratio
    #3. low coverage positions: show coverage

    lowCut = 10 # low cverage for reporting low coverage position
    heteroCut_upper = 0.1
    heteroCut_lower = 0.00001#cut for sequencing errors
    indelCut = 0.1
    types = list('ATCG')
    rdict = {a: {b:0 for b in types} for a in types}
    freq_dict = {a: 0 for a in types}
    comp_count = {a: 0 for a in types}
    with open(fn) as fp:
        for line in fp:
            items = line.strip().split()
            #LEP28439_wenlin_mito    1       A       565     ^].
            try:
                scaf, index, basechar, cov, label, Qlabel = items
            except:
                print('lowCov: %s %s %s %s' % (scaf, index, basechar, cov))
                pass

            if basechar != 'N':
                comp_count[basechar] += 1

            if int(cov) <= lowCut:
                print('lowCov: %s %s %s %s' % (scaf, index, basechar, cov))
                pass

            string = []
            indel = []
            isSkip = False

            label = label.upper().replace('$', '').replace('*', '')
            #indel_flag = 0
            each_indel = []
            pointer = 0
            while(pointer < len(label)):
                char = label[pointer]
                pointer += 1
                if char == '^':
                    pointer += 1
                    continue
                if char == '*':
                    continue
                if char == '-' or char == '+':
                    each_indel = [char]
                    indel_numb = []
                    while(label[pointer].isdigit()):
                        indel_numb.append(label[pointer])
                        pointer += 1

                    each_indel += indel_numb
                    indel_numb = int(''.join(indel_numb))
                    each_indel.append(label[pointer:pointer+indel_numb])
                    pointer += indel_numb
                    indel.append(''.join(each_indel))
                    continue

                string.append(char)

            string = ''.join(string)
            check = string.replace(',', '').replace('.', '')
            if len(check) > (heteroCut_upper * len(string)):
                #not low frequency positions, continue
                continue

            else:
                freq_dict[basechar] += string.count(',') + string.count('.')
                count_dict = Counter(check)
                goodChars = set([char for char in count_dict
                    if count_dict[char] > 2 and (count_dict[char] >= (len(string) * heteroCut_lower))])
                for char in check:
                    if char in goodChars:
                        rdict[basechar][char] += 1

    allN = float(sum(comp_count.values()))
    norm_dict = {char: comp_count[char]/allN*10 for char in comp_count}
    #print 'norm_dict', norm_dict
    header = 'letter\t%s' % ('\t'.join(types))
    print(header)
    for a in types:
        chars1 = [a, reverse_chars[a]]
        line = [a]
        for b in types:
            chars2 = [b, reverse_chars[b]]
            freq = 0
            for char1, char2 in zip(chars1, chars2):
                norm1 = float(freq_dict[char1])
                freq += rdict[char1][char2] / norm1
            line.append(freq)
        line = '\t'.join(map(str, line))
        print(line)
