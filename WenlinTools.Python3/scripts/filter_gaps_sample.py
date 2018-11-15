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
    fn = sys.argv[1]
    sampleIDs = set(cmn.getid(sys.argv[2]))
    print(sampleIDs)

    gapped = set([])
    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
            else:
                seq = line.strip()
                adict[defline] = seq
                if defline.split('_')[0][1:] in sampleIDs:
                    count = 0
                    for char in seq:
                        if char == '-' or char == 'N' or char == ',' or char == 'X':
                            gapped.add(count)
                        count += 1

    mapInfo = ['#both the index started with 0 (not 1)']
    mapDict = {}
    new = []
    for defline in adict:
        seq = adict[defline]
        label = defline[1:].split('_')[0]
        #goodSeq = [char for i, char in enumerate(seq)
        #        if i not in gapped]
        count = 0
        goodSeq = []
        for i, char in enumerate(seq):
            if i not in gapped:
                #mapInfo.append('%s\t%s\t%s' % (label, count, i))
                mapDict[count] = i
                count += 1
                goodSeq.append(char)

        goodSeq = ''.join(goodSeq)

        new.append('%s\n%s\n' % (defline, goodSeq))

    f_label = '%s_%s' % (cmn.lastName(fn).split('.')[0], cmn.lastName(sys.argv[2]))
    dn = f_label + '_noGap.fa'
    cmn.write_file(''.join(new), dn)

    mapInfo += ['%s\t%s' % (key, mapDict[key]) for key in mapDict]
    dn = f_label + '_noGap2coding_index.info'
    cmn.write_lines(mapInfo, dn)

