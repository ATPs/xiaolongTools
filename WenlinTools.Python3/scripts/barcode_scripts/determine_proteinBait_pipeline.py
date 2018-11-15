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
    #try:
    #    fn=sys.argv[1]
    #except:
    #    fn = 'checkSummary.report'
    fn = 'final_comparison.table'


    diffN_dict = {}
    gapDict = {}
    for line in cmn.file2lines(fn)[1:]:
        items = line.strip().split()
        sp = items[0].split('_')[0]
        diffN = int(items[3])
        lv = items[6]
        if 'onfident' in lv:
            continue
        if 'goodRef' in lv:
            continue
        if 'denovo' in lv:
            continue

        seq = items[-1]

        Ngap = seq.count('-') + seq.count('N')
        try:
            gapDict[sp] = min(Ngap, gapDict[sp])
        except:
            gapDict[sp] = Ngap

        if Ngap > 5:
            continue
        try:
            diffN_dict[sp] = min(diffN_dict[sp], diffN)
        except KeyError:
            diffN_dict[sp] = diffN

    badSp = set([sp for sp in diffN_dict if diffN_dict[sp] > 10])
    badSp2 = set([sp for sp in gapDict if gapDict[sp] > 5])
    checkSet = (badSp | badSp2)
    if len(checkSet) > 0:
        sp = os.getcwd().strip('/').split('_')[-1]
        print('launch proteinbait pipeline for %s' % sp)
        cmd = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/run_proteinBait_pipeline.py fqlist %s;' % sp
        cmd += 'bash /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/error_detection.sh %s' % sp
        cmn.run(cmd)
    else:
        print('no need for protein bait for %s' % sp)




