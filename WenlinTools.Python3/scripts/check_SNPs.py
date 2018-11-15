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
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    count = 0
    SNP_count = 0
    heter_count = 0
    with open(fn) as fp:
        for line in fp:
            count += 1
            if 'HaplotypeScore' in line:
                SNP_count += 1
                if '0/1' in line:
                    heter_count += 1

    fraction = float(SNP_count) / count
    f2 = float(heter_count) / count
    print(fn, SNP_count, count, fraction, f2)


