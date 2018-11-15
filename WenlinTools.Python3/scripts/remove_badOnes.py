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
        fn, fadd = sys.argv[1:3]
    except:
        print("Usage: *.py aln.fa repID.file", file=sys.stderr)
        sys.exit()

    #fID = '/work/biophysics/wli/introgression2/4_filterIntro/rep_sps'
    #fID = '/project/biophysics/Nick_lab/wli/sequencing/myAnalysis/clean_ref_bias/4_build_tree/pureIDs'
    #goodIDs = set([i.split()[0] for i in cmn.file2lines(fID)
    #       if i.strip() != ''])

    #fadd = 'added_sps'
    #if cmn.filexist(fadd):
    #    print 'found local list, add them in'
    badIDs = set([each.split()[0].split('_')[0] for each in cmn.file2lines(fadd)
        if each[0] != '#'])

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '_filtered.fa'
    dp = open(dn, 'w')

    new = []
    #leftIDs = set(badIDs)
    leftIDs = set([])
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].strip().strip().split('_')[0].replace('flt', '')
                if name not in badIDs:
                    isGood = True
                else:
                    leftIDs.add(name)
                    isGood = False

            if isGood:
                dp.write(line)

    dp.close()
    if len(leftIDs) != 0:
        print('removed IDs: ')
        print('\n'.join(leftIDs))
    else:
        print('remove nothing!')

