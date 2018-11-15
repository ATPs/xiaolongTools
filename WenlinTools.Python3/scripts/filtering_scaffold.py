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
        lengthCut = int(sys.argv[2])
    except:
        print("Usage: *.py ass.fa 2000", file=sys.stderr)
        sys.exit()


    dn = '%s_cut%s.fa' % (cmn.lastName(fn).replace('.fa', ''), lengthCut)
    dp = open(dn, 'w')

    seq = []
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '>':
                #output once
                if seq != []:
                    seq = ''.join(seq)
                    if len(seq) >= lengthCut:
                        fasta = '>%s\n%s\n' % (name, seq)
                        dp.write(fasta)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
    if seq != []:
        seq = ''.join(seq)
        if len(seq) >= lengthCut:
            fasta = '>%s\n%s\n' % (name, seq)
            dp.write(fasta)

    dp.close()


