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




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        fass = sys.argv[2]
    except:
        print("Usage: *.py log ass.fa", file=sys.stderr)
        sys.exit()

    indel_info = []
    with open(fn) as fp:
        for line in fp:
            if line.startswith('indel'):
                #indel: LEP28439_wenlin_mito 5807 A 85 10 -2AT,-2AT,-2AT,-2AT,-2AT,-2AT,-2AT,-2AT,-2AT,-2AT
                items = line.strip().split()
                index = int(items[2])
                char = items[3]
                labels = items[-1].split(',')
                #now assuming
                #1. every line read needed to be changed
                count_dict = Counter(labels)
                change = max(labels, key=lambda x: count_dict[x])
                indel_info.append((index, char, change))


    seqDict = read_fa(fass)
    name = list(seqDict.keys())[0]
    seq = list(seqDict[name])
    for index, char, change in indel_info:
        i = index - 1
        sign = change[0]
        numb = int(change[1])
        newChars = change[2:]
        if sign == '-':
            for nextI in range(numb):
                seq[i+1+nextI] = 'gap'
        elif sign == '+':
            seq[i] += newChars
        else:
            print('Error: werid sign %s' % sign)

    new = ''.join([each for each in seq
        if each != 'gap'])
    dn = 'fixed.fa'
    fasta = '>%s\n%s\n' % (name, new)
    cmn.write_file(fasta, dn)


