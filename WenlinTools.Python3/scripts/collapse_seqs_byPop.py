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
gapChars = set(['N', '-', '*'])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
transfer_dict = {
        'CT': 'R',
        'AG': 'Y',
        'AT': 'W',
        'CG': 'S',
        'GT': 'M',
        'AC': 'K'

        }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
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


def parse_popDef_newFormat(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        line = line.strip()
        if line == '':
            continue
        if line[0] == '#':
            continue

        items = line.strip().split()
        sp, pop = items[:2]
        try:
            adict[pop].add(sp)
        except KeyError:
            adict[pop] = set([sp])
    return adict


if __name__=='__main__':
    #fn = 'filtered_genomes_coding_noGap.fa'
    fn = sys.argv[1]
    fpoplist = sys.argv[2]

    seqDict = read_fa(fn)
    popDict = parse_popDef_newFormat(fpoplist)

    dp = open('%s_collapsed.fa' % cmn.lastName(fn).replace('.fa', ''), 'w')
    for popname in popDict:
        aset = popDict[popname]
        print('collapsing for %s (%s samples)' % (popname, len(aset)))
        goodNames = [name for name in seqDict
                if name.split('_')[0] in aset]
        seqlist = [seqDict[name] for name in goodNames]

        #collapsing letter
        new = []
        for i in range(len(seqlist[0])):
            chars = [seq[i] for seq in seqlist
                    if seq[i] not in gapChars]
            if len(chars) == 0:
                new.append('-')
                continue
            else:
                count_dict = Counter(chars)
                maxChar = max(count_dict, key=lambda x: count_dict[x])
                if count_dict[maxChar] > len(chars) * 0.8:
                    new.append(maxChar)
                else:
                    new.append('N')
        fasta = '>%s\n%s\n' % (popname, ''.join(new))
        dp.write(fasta)
    dp.close()
