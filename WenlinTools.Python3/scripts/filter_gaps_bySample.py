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
gapChars = set(['X','N','-'])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #fn = 'coding.fasta'
    try:
        fn = sys.argv[1]
        #percentage = float(sys.argv[2]) # 0.3
        sampleID = sys.argv[2]
    except:
        print('*.py *.fa sampleID')
        #print '1.0 would accept any positions'
        sys.exit()

    #1. first read, get gap positions
    takenSeq = None

    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
                if sampleID in defline:
                    isTaken = True
                else:
                    isTaken = False
            else:
                if isTaken:
                    takenSeq = line.strip()
                    print('take the base seq as %s  for %s' % (defline, sampleID))
                    break

    goodP = [i for i in range(len(takenSeq))
            if takenSeq[i] not in gapChars]

    f_label = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '')
    dn = f_label + '_base%s.fa' % sampleID
    dp = open(dn, 'w')

    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
            else:
                goodSeq = ''.join([line[i] for i in goodP])
                dp.write('%s\n%s\n' % (defline, goodSeq))
    dp.close()

