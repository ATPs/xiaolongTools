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
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.sam", file=sys.stderr)
        sys.exit()


    rdict = {}
    scaf_lengths = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '@':
                if 'LN' in line:
                    #@SQ     SN:scaffold25_cov66     LN:268251
                    items = line.strip().split()
                    scaf = items[1][3:]
                    length = int(items[2][3:])
                    scaf_lengths[scaf] = length
                continue

            items = line.strip().split()
            scaf = items[2]
            seq = items[9]
            length = len(seq)
            try:
                rdict[scaf] += length
            except KeyError:
                rdict[scaf] = length


    new = []
    for scaf in scaf_lengths:
        length = scaf_lengths[scaf]
        try:
            cov = rdict[scaf]
        except KeyError:
            cov = 0

        cov = float(cov) / length
        new.append('%s\t%s\t%s\n' % (scaf, length, cov))

    dn = cmn.lastName(fn) + '.scafCov'
    cmn.write_file(''.join(new), dn)

