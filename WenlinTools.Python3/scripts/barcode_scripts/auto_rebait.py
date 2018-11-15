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

def parse_sam_for_stat(fn):
    rdict = {}
    samfile = pysam.AlignmentFile(fn)
    for record in samfile:
        if record.is_unmapped or record.is_secondary:
            continue

        scaf = record.reference_name
        aligns = record.get_aligned_pairs()
        N = len([each for each in aligns
            if None not in each])

        if N >= 30:
            try:
                rdict[scaf] += 1
            except KeyError:
                rdict[scaf] = 1
    return rdict



def extract_same_genus(genus, fall):
    dn = 'genus_for_autoPicking.fa'
    namelist = []
    with open(fall) as fp, open(dn, 'w') as dp:
        for line in fp:
            if '>' in line:
                name = line[1:].strip()
                label = name.split('_')[0]
                if label == genus:
                    isGood = True
                    namelist.append(name)
                else:
                    isGood = False
            if isGood:
                dp.write(line)

    cmd = 'module add bwa; bwa index %s' % dn
    cmn.run(cmd)

    return dn, namelist


def find_genus_info():
    genus = None
    try:
        genus = sys.argv[2]
    except:
        fn = 'restricted_genus.info'
        if cmn.filexist(fn):
            genus = cmn.txt_read(fn).strip()
    return genus

if __name__=='__main__':
    #options=parse_options()
    try:
        fn = sys.argv[1]
    except:
        print("Usage: *.py fqlist genus", file=sys.stderr)
        sys.exit()

    #this code has three types of input:
    #   1. if input is only the fqlist, then require the top one to be much
    #       better than the rest
    #   2. if input contains the genus name, then pick out the genus and then
    #       pick the top one. in this case, don't need the top one to be much
    #       better
    #   3. if input has no genus name but there is a file named
    #       'restricted_genus.info', then the code would take this genus and
    #       then perform things similar to #2

    #decide whether the top hit must be strictly good
    isStrict = False
    fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'

    genus = find_genus_info()

    if genus != None:
        print('found genus constraint as %s, only find barcodes from this genus' % genus)
        fsearch, namelist = extract_same_genus(genus, fall)
        if len(namelist) == 1:
            cmn.write_file(namelist[0], 'picked_bait.txt')
            sys.exit()
    else:
        isStrict = True
        print('use strict criterion for all barcodes')
        fsearch = fall

    statDict = {}
    for fq in cmn.file2lines(fn):
        print('running bwa for %s' % fq)
        label = cmn.lastName(fq)
        fsam = '%s_pickBait.sam' % label
        cmd = 'module add bwa;bwa mem -t 48 -M %s %s > %s 2> %s_bwaTmp.log' % (fsearch, fq, fsam, fsam)
        cmn.run(cmd)
        subDict = parse_sam_for_stat(fsam)
        for key in subDict:
            count = subDict[key]
            try:
                statDict[key] += count
            except KeyError:
                statDict[key] = count


    names = sorted(statDict, key=lambda x: statDict[x], reverse=True)
    lastCount = -1
    takenNames = []
    isTaken = False

    takenNames = []
    if isStrict:
        bestName = names[0]
        if statDict[bestName] > 5 * statDict[names[1]] and (statDict[bestName] - statDict[names[1]] > 100):
            takenNames = [bestName]
    else:
        takenNames = [names[0]]
        bestcount = statDict[names[0]]
        for name in names[1:3]:
            if count >= 0.8 * lastCount:
                takenNames.append(name)

    with open('mapping_stat.info', 'w') as dp:
        for name in names:
            count = statDict[name]
            line = '%s\t%s\n' % (count, name)
            dp.write(line)

    if len(takenNames) == 0:
        cmn.run('rm picked_bait.txt 2> /dev/null')
    else:
        cmn.write_file('\n'.join(takenNames), 'picked_bait.txt')



