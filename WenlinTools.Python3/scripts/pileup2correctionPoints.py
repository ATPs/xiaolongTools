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
    heteroCut = 0.05
    indelCut = 0.1
    with open(fn) as fp:
        for line in fp:
            items = line.strip().split()
            #LEP28439_wenlin_mito    1       A       565     ^].
            try:
                scaf, index, basechar, cov, label, Qlabel = items
            except:
                print('lowCov: %s %s %s %s' % (scaf, index, basechar, cov))
                pass
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
            if len(check) > (heteroCut * len(string)):
                if basechar == 'N':
                    continue

                count_dict = Counter(check.upper())
                maxChar = max(count_dict, key=lambda x: count_dict[x])
                print('hetero: %s %s %s %s %s %s %s' % (scaf, index, basechar, maxChar, cov, float(len(check))/len(string), ','.join(indel)))

            #if len(indel) != 0:
                #print 'indel check', indel
            if len(indel) > (indelCut * len(string)):
                print('indel: %s %s %s %s %s %s' % (scaf, index, basechar, cov, len(indel), ','.join(indel)))
