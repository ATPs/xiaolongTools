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


def sam2mapped_IDs(fsam):
    aset = set([])
    for record in pysam.AlignmentFile(fsam):
        if record.is_unmapped:
            continue
        if record.is_secondary:
            continue
        aset.add(record.query_name)
    return aset



if __name__=='__main__':
    fqlist = 'fqlist'

    fqlist = cmn.file2lines(fqlist)
    if True:
        if len(fqlist) == 3:
            R1 = [each for each in fqlist if '_R1' in each][0]
            R2 = [each for each in fqlist if '_R2' in each][0]
            single = [each for each in fqlist if '_single' in each][0]
        elif len(fqlist) == 2:
            R1 = [each for each in fqlist if '_R1' in each][0]
            R2 = [each for each in fqlist if '_R2' in each][0]
            single = None
        else:
            print('don\'t know how to pair %s' % str(fqlist))

    paired_sam = 'paired.sam'

    mapped_IDs = sam2mapped_IDs(paired_sam)
    print(len(mapped_IDs))

    dn1 = 'forDenovo_R1.fastq'
    dn2 = 'forDenovo_R2.fastq'
    with open(R1) as fp1, open(R2) as fp2, open(dn1, 'w') as dp1, open(dn2, 'w') as dp2:
        for i, line1 in enumerate(fp1):
            line2 = fp2.readline()
            if i % 4 == 0:
                ID = line1[1:].split()[0].split('/')[0]
                if ID in mapped_IDs:
                    isGood = True
                else:
                    isGood = False
            if isGood:
                dp1.write(line1)
                dp2.write(line2)

    if single == None:
        cmn.run('touch forDenovo_singleton.fa')
    else:
        with open(single) as fp, open('forDenovo_singleton.fa', 'w') as dp:
            for i, line in enumerate(fp):
                if i % 4 == 0:
                    ID = line1[1:].split()[0].split('/')[0]
                    if ID in mapped_IDs:
                        isGood = True
                        dp.write('>' + line[1:].replace(' ', '_'))
                    else:
                        isGood = False

                elif i % 4 == 2:
                    if isGood:
                        dp.write(line)
