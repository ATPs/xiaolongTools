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
        print("Usage: *.py vcf", file=sys.stderr)
        sys.exit()


    total = cmn.cmd2info('wc -l %s' % fn).split()[0]

    SNPs = cmn.cmd2info('grep HaplotypeScore %s > %s.tmp; wc -l %s.tmp' % (fn, fn, fn)).split()[0]

    lowqual = cmn.cmd2info('grep LowQual %s.tmp|wc -l ; rm %s.tmp' % (fn, fn)).split()[0]

    print(cmn.lastName(fn), total, SNPs, lowqual, int(SNPs)/float(total), int(lowqual)/float(SNPs))



