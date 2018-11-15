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
        fns = sys.argv[1:]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    fqs = []
    for fn in fns:
        fqs += cmn.getid(fn)


    sequences = []
    for fq in fqs:
        print('reading %s' % fq)
        with open(fq) as fp:
            for i, line in enumerate(fp):
                if i % 4 == 1:
                    sequences.append(line.strip())


    maxlength = max(list(map(len, sequences)))

    for i in range(3, maxlength):
        print('checking %s' % i)
        sample_times = 10
        check_seeds = generate_random_sequence(i, sample_times)
