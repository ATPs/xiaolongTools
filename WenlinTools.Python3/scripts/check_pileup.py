#!/usr/bin/env python

#function: check consistency in pileup file
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import re
from collections import Counter

def remove_indel(chars):
    new = []
    i = 0
    length = len(chars)
    while(i < length):
        char = chars[i]
        if char != '+' and char != '-':
            new.append(char)
        else:#find indel
            #find how many digit
            i += 1
            numb = []
            while(chars[i].isdigit()):
                numb += chars[i]
                i += 1
            numb = int(''.join(numb))
            i += numb - 1

        i += 1
    return ''.join(new)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_read_stack(stack, ref_base, count):
    #removed_chars = '. , ^] $ ^O A T G C *'.split()
    #re_chars = ['\^.{1}', '-[0-9]+[ACGTNacgtn]+', '\+[0-9]+[ACGTNacgtn]+', '\$',]
    re_chars = ['\^.{1}', '\$',]
    check_chars = stack.upper()

    for char in re_chars:
        check_chars = re.sub(char, '', check_chars)

    check_chars = remove_indel(check_chars)

    #if check_chars.strip() != '':
    #    print count, stack
    #    print check_chars
    #    sys.exit()

    # . , ATCG should be total count
    if len(check_chars) != int(count):
        print(count, len(check_chars), check_chars)
        print(line)
        sys.exit()

    check_chars = check_chars.replace(',', '.')
    check_chars = check_chars.replace('.', ref_base)

    countDict = Counter(check_chars)
    #stat = '|'.join(['%s:%s' % (i, countDict[i]) for i in countDict])
    return countDict


def count2stat(adict):
    alist = [('%s:%s' % (i, adict[i]), adict[i])
            for i in adict]

    slist = sorted(alist, cmp=lambda x, y: cmp(x[-1], y[-1]), reverse=True)
    stat = '|'.join([i[0] for i in slist])
    return stat


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py 1708_mapped.pileup", file=sys.stderr)
        sys.exit()

    new = []
    try:
        coding = cmn.pickle_read('coding.indexes.pkl')
    except:
        print('do not find index file for coding region')
        print('would not label coding positions')
        coding = set([])


    with open(fn) as fp:
        for line in fp:
            try:
                scaffold, index, ref_base, count, read_stack, qual_stack  = line.strip().split()
            except:
                #coverage is 0
                items = line.strip().split()
                index, ref_base, count = items[1:4]
                print(index, ref_base, count, '0 lowCoverage(0)')
                continue
            countDict = parse_read_stack(read_stack, ref_base, count)

            label = ''
            if int(index) in coding:
                label += 'coding,'

            count = int(count)
            if '*' in countDict:
                count -= countDict['*']
                label += 'gap:%s,' % countDict['*']
            del countDict['*']

            ref_count = countDict[ref_base]
            #if ref_count < count * 0.6:
            if count < 10:
                label += 'lowCoverage,'
            if ref_count < count * 0.7:
                label += 'badP,'

            #check polymophism
            polyM = []
            for char in countDict:
                char_count = countDict[char]
                if char_count > 50:# need to have enough coverage
                    if char_count > (0.1 * count):#more than 10%
                        polyM.append(char)
            if len(polyM) > 1:
                label += 'polyM:%s,' % ('-'.join(polyM))

            if count == 0:
                freq = 0
            else:
                freq = float(countDict[ref_base]) / count
                freq = '%.2f' % freq
            stat = count2stat(countDict)
            print(index, ref_base, count, freq, stat, label)




