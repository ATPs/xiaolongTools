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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #fn = 'coding.fasta'
    try:
        fn = sys.argv[1]
        percentage = float(sys.argv[2]) # 0.3
    except:
        print('*.py *.fa 0.3')
        print('1.0 would accept any positions')
        sys.exit()

    #1. first read, get gap positions
    gapped = {}
    seqCount = 0
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
                seqCount += 1
            else:
                count = 0
                seq = line.strip()
                for char in seq:
                    if char == '-' or char == 'X' or char == '*':
                        try:
                            gapped[count] += 1
                        except:
                            gapped[count] = 1
                    count += 1

    length = len(seq)
    #2. tell which positions is not gapped
    cutoff = percentage * seqCount
    badP = set([i for i in gapped if gapped[i] > cutoff])
    #print goodP

    #mapInfo = ['#both the index started with 0 (not 1)']
    #mapDict = {}
    new = []
    #total = len(adict)
    f_label = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '')
    dn = f_label + '_%sGap.fa' % percentage
    dp = open(dn, 'w')

    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
            else:
                goodSeq = ''.join([line[i] for i in range(length) if i not in badP])
                dp.write('%s\n%s\n' % (defline, goodSeq))
    dp.close()

