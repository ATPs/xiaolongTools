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
import re


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    order = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0].strip()
        seq = ''.join(lines[1:])
        seq = list(seq.strip())
        #Ngap = seq.count('-')
        #if Ngap > (0.8 * len(seq)):
        #    continue
        adict[defline] = seq
        order.append(defline)
    return adict, order

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def scaf2numb(scaf):
    numb = cmn.find_between(scaf, 'scaffold', '_cov')
    return int(numb)



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fr = sys.argv[1:3]
    except:
        print("Usage: *.py realigned.pileup assembly_v2.fa", file=sys.stderr)
        sys.exit()


    adict = {}

    bp = 'A T G C'.split()

    with open(fn) as fp:
        for line in fp:
            #scaffold1_cov51 1       A       2       ^].^].  II
            items = line.strip().split()
            scaf = items[0]
            index = int(items[1])
            #the pileup file would miss the positions with no reads
            if scaf not in adict:
                adict[scaf] = []
                counter = 1

            while(index != counter):
                counter += 1
                adict[scaf].append('N')

            #now the index meets the counter, tell the dominated form
            base = items[2]
            total_count = int(items[3])

            if total_count == 0:
                counter += 1
                adict[scaf].append('N')
                continue

            stat = items[4]

            stat = stat.replace(',', '.')
            stat = stat.replace('.', base)
            stat = re.sub('-[0-9]+[ACGTNacgtn]+', '', stat)
            stat = re.sub('\+[0-9]+[ACGTNacgtn]+', '', stat)

            stat = stat.upper()

            maxBase = (base, 0)
            for char in bp:
                count = stat.count(char)
                if count > maxBase[1]:
                    maxBase = (char, count)
                else:
                    #here is equal or less
                    if count == maxBase[1] and char != base:
                        maxBase = (char, count)

            adict[scaf].append(maxBase[0])
            counter += 1

    seqDict, scafs = read_fa(fr)

    #output the dict
    fastas = []

    #scafs = seqDict.keys()
    #scafs = sorted(scafs, cmp=(lambda x, y: cmp(scaf2numb(x), scaf2numb(y))))

    for scaf in scafs:
        try:
            seq1 = adict[scaf]
        except:
            seq1 = []

        seq2 = seqDict[scaf]

        if len(seq1) != len(seq2):
            print('warning: the length of %s has been changed!' % scaf)

        newSeq = []
        for i, char2 in enumerate(seq2):
            try:
                char1 = seq1[i]
                if char1 != 'N' and char1 != '-':
                    newSeq.append(char1)
                else:
                    newSeq.append(char2)
            except:
                newSeq.append(char2)

        fastas.append('>%s\n%s\n' % (scaf, ''.join(newSeq)))
    dn = 'assembly_selfref_v2.fa'
    cmn.write_file(''.join(fastas), dn)
