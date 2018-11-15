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
        print("Usage: *.py fastq", file=sys.stderr)
        sys.exit()

    Min = 'h'
    Max = '!'

    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 3:
                qletters = line.strip()
                Min = min(min(qletters), Min)
                Max = max(max(qletters), Max)

    print('Max Value:', Max, ord(Max))
    print('Min Value:', Min, ord(Min))


