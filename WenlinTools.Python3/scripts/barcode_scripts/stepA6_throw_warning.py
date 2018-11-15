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
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        wdir=sys.argv[1]
    except:
        print("Usage: *.py wdir", file=sys.stderr)
        sys.exit()

    os.chdir(wdir)


    fns = cmn.cmd2lines('ls sampleRun_*/final_comparison.table')

    infoDict = {}
    for line in cmn.file2lines('compare.check'):
        sp = line.strip().split()[0]
        infoDict[sp] = line

    sum_header = None
    sum_content = []

    for fn in fns:
        lines = cmn.file2lines(fn)
        sum_header = lines[0]
        sum_content += lines[1:]


    hasBad = False
    with open('checkSummary.report', 'w') as dp, open('checkWarnning.list', 'w') as dp_warn, open('wenlinWarnning.list', 'w') as dpWenlin:
        dp.write(sum_header + '\n')
        for line in sum_content:
            dp.write(line + '\n')

            items = line.strip().split()
            annotated = items[1].split('_')[0].lower()
            if annotated == 'j.':
                annotated = 'junonia'
            barcodeSP = items[2]
            diffN = float(items[3])
            label = items[0]

            if annotated not in barcodeSP.lower() and diffN < 10:
                if 'Suspicius' in line:
                    hasBad = True
                    dp_warn.write('possibleMessup\t%s\t%s->"%s"\n' % (label.replace('_','\t'), items[1], barcodeSP))
                    dpWenlin.write('possibleMessup\t%s\t%s\n' % (label.replace('_','\t'), infoDict[label.split('_')[0]]))
                elif 'onfident' in line:
                    hasBad = True
                    dp_warn.write('sureMessup\t%s\t%s->%s\n' % (label.replace('_','\t'), items[1], barcodeSP))
                    dpWenlin.write('sureMessup\t%s\t%s\n' % (label.replace('_','\t'), infoDict[label.split('_')[0]]))
                elif 'TODO' in line:
                    seq = items[-1]
                    if seq.count('N') > 0.5 * len(seq):
                        #poor denovo
                        continue
                    dp_warn.write('TODO_check\t%s\t%s->%s\n' % (label.replace('_','\t'), items[1], barcodeSP))
                    dpWenlin.write('TODO: %s\t%s\n' % (label.replace('_','\t'), infoDict[label.split('_')[0]]))


    if hasBad:
        print('something is wrong, please send checkWarnning.list and checkSummary.report to Wenlin and Nick')
        print('remember to rename these two files according to the batch name')
    else:
        print('so far so good')





