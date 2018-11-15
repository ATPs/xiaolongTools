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
import random

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])

        adict[defline] = seq
    return adict, len(seq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':

    #fn = 'all_genomes_noGap.fa'
    #fn = 'all_genomes_charGap.fa'
    try:
        fn, parts = sys.argv[1:3]
        parts = int(parts)
    except:
        print('*.py all_genomes_charGap.fa 100')
        print('this is random split with shuffling')
        sys.exit()

    adict, length = read_fa(fn)

    #random sample
    indexes = list(range(length))
    random.shuffle(indexes)

    #parts = 100
    each = length / parts
    print('randomly split into %s parts (%s bp in each part)' % (parts, each))

    fnlabel = cmn.lastName(fn).replace('.fa', '')

    outdir = 'split_%s' % fnlabel
    cmn.mkdir(outdir)

    keys = list(adict.keys())

    for i in range(0, length, each):
        sub_indexes = indexes[i: i+each]
        subset = []

        for key in adict:
            alist = adict[key]
            seq = ''.join([alist[i] for i in sub_indexes])
            fasta = '>%s\n%s\n' % (key, seq)
            subset.append(fasta)

        dn = '%s/%s_%s.fa' % (outdir, fnlabel, i)
        cmn.write_file(''.join(subset), dn)

