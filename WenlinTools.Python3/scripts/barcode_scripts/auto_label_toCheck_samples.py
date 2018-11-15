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
def find_sp(line):
    sp = line[1:].split('_')[0]
    if ']' in sp:
        sp = sp.split(']')[-1]
    sp = sp.replace('=', '')
    return sp


def old_get_goodIDs():
    fn1 = 'compare.check'
    #fn2 = 'checkSummary.report'
    goodIDs = set([])
    for line in cmn.file2lines(fn1):
        items = line.strip().split()
        if 'takenD' in line:
            goodIDs.add(items[0])

    return goodIDs

def get_goodIDs():
    fn1 = 'compare.check'
    adict = {}
    for line in cmn.file2lines(fn1):
        items = line.strip().split()
        if 'takenD' in line:
            continue
        sp = items[0]
        if 'completeD' in line and ('takenD' not in line):
            adict[sp] = 'likelyWrong'
        elif 'goodCC' in line:
            adict[sp] = 'mightOK'
        else:
            adict[sp] = 'notSure'
    return adict




if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    dn = 'labeled.fa'

    #badIDs = set(cmn.file2lines('badIDs'))
    goodIDs = get_goodIDs()

    taken = set([])
    with open(fn) as fp, open(dn, 'w') as dp:
        for line in fp:
            if '(assembled)' not in line:
                dp.write(line)
            else:
                sp = find_sp(line)
                try:
                    label = goodIDs[sp]
                    if True:
                        print('label %s as %s' % (sp, label))
                        line = '>[%s]%s' % (label, line[1:])
                    taken.add(sp)
                except KeyError:
                    pass
                dp.write(line)

    #missing = badIDs - taken

    #if len(missing) != 0:
    #    print 'These badIDs are not labeled in the file:'
    #    print '\n'.join(missing)


