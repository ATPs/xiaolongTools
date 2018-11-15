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
def scaf2numb(scaf):
    numb = cmn.find_between(scaf, 'scaffold', '_cov')
    return int(numb)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py realigned.pileup", file=sys.stderr)
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
                else:#equal or less
                    if count == maxBase[1] and char != base:#try to be as different as base
                        maxBase = (char, count)

            adict[scaf].append(maxBase[0])
            counter += 1


    #output the dict
    scafs = list(adict.keys())
    scafs = sorted(scafs, cmp=(lambda x, y: cmp(scaf2numb(x), scaf2numb(y))))
    fastas = ['>%s\n%s\n' % (scaf, ''.join(adict[scaf]))
            for scaf in scafs]
    dn = 'assembly_selfref.fa'
    cmn.write_file(''.join(fastas), dn)
