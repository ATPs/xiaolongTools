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

    gapped = {}
    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
            else:
                seq = line.strip()
                adict[defline] = seq
                count = 0
                for char in seq:
                    if char == '-' or char == 'N' or char == ',' or char == 'X':
                        try:
                            gapped[count] += 1
                        except:
                            gapped[count] = 1
                    count += 1

    mapInfo = ['#both the index started with 0 (not 1)']
    mapDict = {}
    new = []
    total = len(adict)
    cutoff = percentage * total
    for defline in adict:
        seq = adict[defline]
        label = defline[1:].split('_')[0]
        #goodSeq = [char for i, char in enumerate(seq)
        #        if i not in gapped]
        count = 0
        goodSeq = []
        for i, char in enumerate(seq):
            try:
                Ngap = gapped[i]
                if Ngap <= cutoff:
                    isGood = True
                else:
                    isGood = False
            except:#no gap
                isGood = True

            if isGood:
                #mapInfo.append('%s\t%s\t%s' % (label, count, i))
                mapDict[count] = i
                count += 1
                goodSeq.append(char)

        goodSeq = ''.join(goodSeq)

        new.append('%s\n%s\n' % (defline, goodSeq))

    f_label = cmn.lastName(fn).split('.')[0]
    dn = f_label + '_%sGap.fa' % percentage
    cmn.write_file(''.join(new), dn)

    mapInfo += ['%s\t%s' % (key, mapDict[key]) for key in mapDict]
    dn = f_label + '_%sGap2coding_index.info' % percentage
    cmn.write_lines(mapInfo, dn)

