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
from random import randint


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_length_info(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        scaf, length = line.strip().split()
        adict[scaf] = int(length)
    return adict


def fill_in_missing(smallI, largeI):
    #fill till the index is the same (the same one is not included)
    while(smallI != largeI):
        print('N\tN')
        smallI += 1

def append_last_chars(smallI, length):
    if smallI != length:
        #if the length is not meet, fill to the length
        for _ in range(length - smallI):
            print('N\tN')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def vcf2map(line):
    items = line.strip().split()
    scaf, index = items[:2]
    base, mut = items[3:5]

    info = items[-1]
    subitems = info.split(':')
    GT = subitems[0]
    #if 'LowQual' in line:
        #not sufficient evident to call anything
        #keep original ones
    #    fasta[0].append(base)
    #    fasta[1].append(base)
    #    continue

    char1, char2 = '', ''
    if GT == './.':#no coverage
        if items[7] == '.':
            char1 = '-'
            char2 = '-'
        else:
            char1 = 'N'
            char2 = 'N'
        #print '%s\t%s' % (char1, char2)
        #return [char1, char2]
        return '%s\t%s' % (char1, char2)

    GT_counts = subitems[1].split(',')

    if GT == '0/0': #no SNP call
        char1 = base
        char2 = base
    elif GT == '1/1': #totally different
        char1 = mut
        char2 = mut

    elif GT == '0/1': #called SNP
        i, j = list(map(int, GT_counts))
        randomN = randint(0, 1)
        char1 = base
        char2 = mut
        if randomN:
            char1, char2 = char2, char1

    elif GT == '1/2':
        o, i, j = list(map(int, GT_counts))
        if o != 0:
            print('warning: original frequency is not 0 for:', file=sys.stderr)
            print(line, file=sys.stderr)
        char1, char2 = mut.split(',')
        randomN = randint(0, 1)
        if randomN:
            char1, char2 = char2, char1

    else:
        print('unrecognized line: %s' % line, file=sys.stderr)

    #return [char1, char2]
    return '%s\t%s' % (char1, char2)



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fdict = sys.argv[1:3]
    except:
        print("Usage: *.py *.vcf length.txt", file=sys.stderr)
        sys.exit()

    scaf_length = read_length_info(fdict)

    #used to record the outputted index
    scaf_counters = {}

    #used to make sure the last scaffold output to the end
    lastIndex = -1
    lastScaf = ''

    with open(fn) as fp:
        for line in fp:
            #if not line.startswith('scaffold'):
            if line.strip()[0] == '#':
                continue

            items = line.strip().split()
            scaf = items[0]

            try:
                expect_index = scaf_counters[scaf]
            except:
                expect_index = 1
                #a new scaffold, check if it output to the end
                try:
                    append_last_chars(lastIndex, scaf_length[lastScaf])
                except:
                    if lastScaf != '':
                        print('warning: scaffold %s not found!' % lastScaf, file=sys.stder)

            index = int(items[1])
            #make the expect index run the the current
            fill_in_missing(expect_index, index)

            #parse current
            print(vcf2map(line))

            scaf_counters[scaf] = index + 1
            lastScaf = scaf
            lastIndex = index


    #append for the last scaffold in the list
    append_last_chars(lastIndex, scaf_length[lastScaf])




