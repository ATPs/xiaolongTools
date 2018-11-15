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
from fractions import Fraction


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_percent(percent):
    numb = Fraction(percent)
    a = numb.denominator
    b = numb.numerator
    return a, b


if __name__=='__main__':
    #options=parse_options()
    try:
        fR1, fR2 = sys.argv[1:3]
        percent = sys.argv[3]
    except:
        print("Usage: *.py R1 R2 0.8", file=sys.stderr)
        print('0.8 (80%) is the percentage left after filtering', file=sys.stderr)
        sys.exit()


    sample = cmn.lastName(fR1).split('_')[0]

    denominator, numerator = parse_percent(percent)
    numerator -= 1

    count = 0
    fpR1 = open(fR1)
    fpR2 = open(fR2)
    dnlabel = '%s_p%s' % (sample, percent)
    print('making %s' % dnlabel)

    dnR1 = open('%s_R1.fastq' % dnlabel, 'w')
    dnR2 = open('%s_R2.fastq' % dnlabel, 'w')

    for i, line1 in enumerate(fpR1):
        line2 = fpR2.readline()
        if i % 4 == 0:
            count += 1

        if count % denominator <= numerator:
            dnR1.write(line1)
            dnR2.write(line2)

    dnR1.close()
    dnR2.close()
    fpR1.close()
    fpR2.close()



