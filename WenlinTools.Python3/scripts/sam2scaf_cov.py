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
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py sam", file=sys.stderr)
        sys.exit()


    samfile = pysam.AlignmentFile(fn)

    len_dict = {}
    for record in samfile.header['SQ']:
        len_dict[record['SN']] = record['LN']

    count_dict = {}
    for read in samfile:
        if read.is_unmapped or read.is_secondary:
            continue

        scaf = read.reference_name
        for i,j in read.aligned_pairs:
            if j != None:
                try:
                    count_dict[scaf] += 1
                except KeyError:
                    count_dict[scaf] = 1


    new = []
    for scaf in count_dict:
        count = count_dict[scaf]
        try:
            length = len_dict[scaf]
        except:
            print('warning: can not find %s in header, skip' % scaf)

        line = [scaf, length, float(count) / length]
        new.append('\t'.join(map(str, line)) + '\n')

    dn = cmn.lastName(fn).replace('.sam', '') + '_cov.txt'
    cmn.write_file(''.join(new), dn)




