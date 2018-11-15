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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main

def isHighRatio(line):
    ratio = float(line.strip().split()[5].split('(')[0])
    if ratio > 0.3:
        return True
    else:
        return False

def read_compare_file(fn, fcheck, ID):
    cmd = 'grep -P \'^%s\t\' %s' % (ID, fn)
    line = cmn.cmd2lines(cmd)[0]
    #TODO: need to rework
    if 'same' in line:
        return 'mostConfident', [0]
    elif 'diffGenus' in line:
        return 'diffGenus', [0, 1]
    elif 'takenD' in line:
        return 'confident', [0]
    elif 'completeDenovo' in line and ('goodCC' in line):
        return 'confident', [0]
    elif 'completeDenovo' in line:
        return 'denovoOnly', [0]
    elif 'goodCC' in line:
        return 'goodRef', [1]
    else:
        #these sample would be somewhat problematic
        if 'Error' in line:
            return 'Error', [0, 1]
        elif isHighRatio(line):
            return 'Suspicius', [0, 1]
        elif 'Gap0' not in line:
            return 'PoorSample', [1]
        elif 'noDenovo' in line:
            return 'refOnly', [0]
        else:
            return 'TODO', [0, 1]


def isConsistent(seq1, seq2):
    gapChars = set(['N', '-'])
    for char1, char2 in zip(seq1, seq2):
        if char1 in gapChars:
            continue
        if char2 in gapChars:
            continue
        if char1 != char2:
            return False
    return True

def check_prot_consistency(seq0, keylist, barcodes):
    flag_list = []
    for key in keylist:
        seqlist = barcodes[key]
        name, seq = seqlist[0]
        label = name.split('_')[-1]
        if isConsistent(seq0, seq):
            flag_list.append('%s_consistent' % label)

    if len(flag_list) == 0:
        return 'unique'
    elif len(flag_list) == len(keylist):
        return 'consistent'
    else:
        return ','.join(flag_list)


def isSameGenus(name1, name2):
    genus = name1.split('_')[0]
    if genus in name2:
        return True
    else:
        return False

if __name__=='__main__':
    sampleID = sys.argv[1]

    cmd = 'cd ..; python /archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/step8_summarize_barcodes.py'
    #if (not cmn.filexist('../compare.check')) or (sampleID not in [each.split()[0] for each in cmn.file2lines('../compare.check')]):
    cmn.run(cmd)

    fexact = '%s.exactCheck' % sampleID
    fcheck = '%s_summaryCheck.table' % sampleID
    fbarcode = '%s_exactCheck.fasta' % sampleID

    confident_level, takenIndex = read_compare_file('../compare.check', fcheck, sampleID)

    dn = 'final_comparison.table'
    dp = open(dn, 'w')
    header = 'orignalID annotatedSP barcodeSP diffN CSid_diffN CSidFrom confidentLevel barcodeSeq'
    dp.write('\t'.join(header.split()) + '\n')

    CSmatch = {}
    for line in cmn.file2lines(fexact):
        sp, fromID, N = line.strip().split('\t')
        #sp = sp.split('_')[0]
        if fromID.strip() == '':
            fromID = 'NA'
        CSmatch[sp] = [fromID, N]


    #sample	checkFlag	diffN	barcodeSP	annotatedSP
    checkDict = {}
    for line in cmn.file2lines(fcheck)[1:]:
        print(line)
        sample, flag, diffN, barcodeSP, annotatedSP = line.strip().split()
        checkDict[sample] = [flag, diffN, barcodeSP, annotatedSP]


    seqDict = read_fa(fbarcode)
    barcodes = {}
    protKey = None
    for name in seqDict:
        #sp = name.split('_')[0]
        sp = name
        fasta = tuple([name, seqDict[name]])
        try:
            barcodes[sp].append(fasta)
        except KeyError:
            barcodes[sp] = [fasta]

        if '_prot' in sp:
            protKey = sp

    keylist = list(barcodes.keys())
    keylist.sort()
    if protKey != None:
        keylist.remove(protKey)
    print('confident_level', confident_level)
    if len(keylist) == 1:
        takenIndex = set(range(len(keylist)))
    for index, sp in enumerate(keylist):
        if index not in takenIndex:
            continue

        line = [sp]
        flag, diffN, barcodeSP, annotatedSP = checkDict[sp]
        annotatedSP = '_'.join(annotatedSP.split('_')[1:])
        fromID, csN = CSmatch[sp]
        line += [annotatedSP, barcodeSP, diffN, csN, fromID]

        seqlist = barcodes[sp]
        if len(seqlist) == 1:
            name, seq = seqlist[0]
            lv = confident_level
        else:
            lv += '_inconsistent'
            seq = ''
            for name, each in seqlist:
                seq += '[%s]%s;' % (name, each)
        #if isSameGenus(annotatedSP, barcodeSP):
        #    lv = 'TODO'
        line += [lv, seq]

        dp.write('\t'.join(line) + '\n')

    if protKey != None:
        sp = protKey
        line = [protKey]
        flag, diffN, barcodeSP, annotatedSP = checkDict[sp]
        annotatedSP = '_'.join(annotatedSP.split('_')[1:])
        fromID, csN = CSmatch[sp]
        line += [annotatedSP, barcodeSP, diffN, csN, fromID]
        seqlist = barcodes[sp]
        name, seq = seqlist[0]
        lv = check_prot_consistency(seq, keylist, barcodes)
        line += [lv, seq]

        dp.write('\t'.join(line) + '\n')

    dp.close()
