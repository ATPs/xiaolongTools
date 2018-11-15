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
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    dn = cmn.lastName(fn) + '.fa'
    samfile = pysam.AlignmentFile(fn)

    dp = open(dn, 'w')

    for record in samfile:
        name = record.query_name
        if record.is_read1:
            name += '_1'
        elif record.is_read2:
            name += '_2'

        seq = record.query_sequence
        fasta = '>%s\n%s\n' % (name, seq)
        dp.write(fasta)

    dp.close()

    print('results are in %s' % dn)



