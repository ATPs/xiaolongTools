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
        print('1.0 would accept any sequence')
        sys.exit()

    new = []
    total = 0
    removeN = 0
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                defline = line.strip()
                total += 1
            else:
                seq = line.strip()
                checkseq = seq.replace('N', '-')
                if checkseq.count('-') <= (percentage * len(seq)):
                    new.append('%s\n%s\n' % (defline, seq))
                else:
                    removeN += 1


    f_label = cmn.lastName(fn).split('.')[0]
    dn = f_label + '_%sGapSample.fa' % percentage
    cmn.write_file(''.join(new), dn)

    print('removed %s from %s sequences' % (removeN, total))
