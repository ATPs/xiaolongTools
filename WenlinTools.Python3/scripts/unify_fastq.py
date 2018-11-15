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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fq", file=sys.stderr)
        sys.exit()


    taken = set([])

    dn = cmn.lastName(fn).split('.')[0] + '_unified.fastq'

    isGood = True
    #new = []
    dp = open(dn, 'w')
    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 0:
                ID = line.strip()
                if ID not in taken:
                    isGood = True
                    taken.add(ID)
                else:
                    isGood = False
                    print('duplcated ID: %s' % ID)
            if isGood:
                dp.write(line)

    #cmn.write_file(''.join(new), dn)
    dp.close()




