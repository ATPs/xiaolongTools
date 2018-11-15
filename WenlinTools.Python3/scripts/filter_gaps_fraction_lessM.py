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
                    if char == '-' or char == 'N' or char == ',' or char == 'X':
                        try:
                            gapped[count] += 1
                        except:
                            gapped[count] = 1
                    count += 1

    #2. tell which positions is not gapped
    cutoff = percentage * seqCount
    goodP = [i for i in gapped if gapped[i] <= cutoff]

    #mapInfo = ['#both the index started with 0 (not 1)']
    #mapDict = {}
    new = []
    #total = len(adict)
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
            else:
                goodSeq = ''.join([line[i] for i in goodP])
                new.append('%s\n%s\n' % (defline, goodSeq))

    f_label = cmn.lastName(fn).split('.')[0]
    dn = f_label + '_%sGap_LM.fa' % percentage
    with open(dn, 'w') as fp:
        for line in new:
            fp.write(line)

