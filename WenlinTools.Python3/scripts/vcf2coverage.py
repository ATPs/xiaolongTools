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
        print("Usage: *.py *.vcf", file=sys.stderr)
        sys.exit()

    #to count the total line
    total_count = 0
    #to save the coverage
    cov_list = []

    with open(fn) as fp:
        for line in fp:
            if line[0] == '#':
                continue

            items = line.strip().split()
            total_count += 1
            codes = items[-2].split(':')
            try:
                i = codes.index('DP')
            except ValueError:
                continue

            values = items[-1].split(':')
            cov_list.append(int(values[i]))

    
    #1. percentage of mapping
    mapped = len(cov_list)
    mPert = float(mapped) / total_count
    cov = sum(cov_list)/float(mapped)
    cov_list.sort()
    cov_M = cov_list[len(cov_list)/2]
    info = 'total length: %s\nmapped percent: %s\ncoverage average: %s\ncoverage median: %s\n' % (total_count, mPert, cov, cov_M)
    info += 'report: %s %s %s %s %s' % (cmn.lastName(fn), total_count, mPert, cov, cov_M)
    #cmn.write_file(info, dn)
    print(info)
