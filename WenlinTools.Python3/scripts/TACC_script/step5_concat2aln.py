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
        print("Usage: *.py result.txt", file=sys.stderr)
        sys.exit()

    #concat results into alignment

    adict = {}
    exon_length = {}
    with open(fn) as fp:
        for line in fp:
            exon, sp, seq = line.strip().split()
            sp = sp.split('.')[0]
            try:
                if exon_length[exon] != len(seq):
                    print('different length found for %s' % line)
            except:
                pass
            exon_length[exon] = len(seq)

            if sp not in adict:
                adict[sp] = {}

            adict[sp][exon] = seq


    new = []
    exons = list(exon_length.keys())
    for sp in adict:
        subdict = adict[sp]
        seq = []
        for exon in exons:
            try:
                segment = subdict[exon]
            except:
                segment = '-' * exon_length[exon]
            seq.append(segment)

        seq = ''.join(seq)
        new.append('>%s\n%s\n' % (sp, seq))

    dn = 'concat_alignment.fa'
    cmn.write_file(''.join(new), dn)



