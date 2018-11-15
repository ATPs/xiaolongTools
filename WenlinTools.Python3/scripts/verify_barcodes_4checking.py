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
from fullname_lib import get_SRNPnumber
import re

SRNP_dict = get_SRNPnumber()


def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0].split()[0].replace('(assembled)', '').strip('.')
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict

def read_fa_barcode(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    alist = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0].replace('flt', '')
        #seq = ''.join(lines[1:]).rstrip('-')
        seq = ''.join(lines[1:])
        if len(seq) > 658:
            seq = seq[20:678]
            #if len(seq) != 698:
            #    print 'warining: length is weird %s (%s)' % (defline, len(seq))
        adict[defline] = seq
        alist.append(defline)
    return adict, alist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def do_barcode_blast(sequence):
    fdb = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_4verify.fa'
    namelabel = sequence.split()[0][1:].split()[0].split('|')[0].replace('*','').split('[')[0].replace('"','').replace("'", '').split('(')[0].split('/')[0].split('\\')[0]
    fquery = '/tmp/%s.fa' % namelabel
    cmn.write_file(sequence, fquery)

    fbr = fquery + '.br'
    cmd = 'module add blast; blastn -max_target_seqs 1000 -query %s -db %s -ungapped ' % (fquery, fdb)
    cmd += '-outfmt \'6 sseqid slen length pident qstart qend qseq sseq\''
    cmd += ' -out %s ' % fbr
    #print cmd
    cmn.run(cmd)
    #cmd += ' | head -n 10'
    #lines = cmn.cmd2lines(cmd)
    lines = cmn.file2lines(fbr)
    cmn.run('rm %s' % fquery)
    return lines

def compare_top_hit(lines):
    items = lines[0].split()
    qstart, qend, qseq, sseq = items[-4:]
    sseqid = items[0]
    qstart = int(qstart)
    adict = {}
    for i, char1 in enumerate(qseq):
        index = i + qstart
        char2 = sseq[i]
        if char1 != char2:
            adict[index] = '%s:%s' % (char1, char2)
    return adict, sseqid


def isGood_ID(ID, name):
    spliter = '[\s|\-_]'
    name = name.replace('-SRNP-', '')
    #print name, spliter
    items = re.split(spliter, name)
    #print ID, items
    if ID in items:
        return True
    return False

def find_barcode_of_itself(ID, br_results):
    global barcode_deflines, barCodeDict
    try:
        SRNP = SRNP_dict[ID]
        print('found SRNP number for %s as %s' % (ID, SRNP))
    except KeyError:
        SRNP = 'somethingNever'

    found_SRNP = False
    noAddback = True
    #print ID, SRNP
    #if the self barcode is inside the list of blast result, make self equal the blast
    for line in br_results:
        items = line.strip().split()
        sseqid = items[0]
        checkList = sseqid.replace('|', '  ').split()
        #print sseqid
        if SRNP in checkList:
            found_SRNP = True

        found_sample = isGood_ID(ID, sseqid)

        if found_sample or found_SRNP:
            sseq = items[-1]
            qstart = int(items[-4])
            qend = int(items[-3])
            if qstart != 1:
                diff = qstart - 1
                sseq = 'N' * diff + sseq

            if qend != 658:
                diff = 658 - qend
                sseq = sseq + 'N' * diff
            #fasta = '>%s\n%s\n' % (sseqid, sseq)
            print('found barcode of itself: %s' % sseqid)
            if sseq.count('N') + sseq.count('-') > 10:#too gapped
                noAddback = False
            return (sseqid, sseq), noAddback

    if not found_SRNP and SRNP != 'somethingNever':
        #try to find by the purename
        Name = [each for each in barcode_deflines
                if SRNP in each.replace('|', '  ').split()]
        if len(Name) != 0:#found it
            defline = Name[0]
            seq = barCodeDict[defline]
            #fasta = '>%s\n%s\n' % (defline, seq)
            print('%s has %s, It is in the verification set but didn\'t find by BLAST..., So add both BLAST best and the exact SRNP number one' % (ID, SRNP))
            return (defline, seq), False

        Name = [each for each in barcode_deflines
                #if any([ID in each.split(separator)
                #    for separator in '|'])]
                if ( ID in each.split('|') or
                     #ID in each.split('_') or
                     ID == each.split('-')[0]
                    )]
        if len(Name) != 0:#found it
            defline = Name[0]
            seq = barCodeDict[defline]
            #fasta = '>%s\n%s\n' % (defline, seq)
            print('found sequence of itself %s for ID %s, It is in the verification set but didn\'t find by BLAST..., So add both BLAST best and the exact SRNP number one' % (defline, ID))
            return (defline, seq), False


        else:
            print('%s has %s, but not found in the verification set' % (ID, SRNP))

    #if not found_sample:
        #try to find the barcode sequence of the sample by the sample ID
        #would only reach here because PCR result generate very partail result and not detectable by blast against database

    return None, False

def fasta2seq(fasta):
    lines = fasta.strip().split('\n')
    seq = ''.join(lines[1:])
    return seq


def compare_itself(fasta1, fasta2):
    seq1 = fasta2seq(fasta1).upper()
    seq2 = fasta2seq(fasta2).upper()
    #seq2 is the barcode
    adict = {}
    index = 0
    gapChar = set(list('-N'))
    for char1, char2 in zip(seq1, seq2):
        index += 1
        if char2 in gapChar:
            continue
        if char1 != char2:
            adict[index] = '%s:%s' % (char1, char2)

    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def rename_pair(name):
    return name.replace('(assembled)', '')


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    seqDict, orderlist = read_fa_barcode(fn)
    barCodeDict = read_fa('/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_4verify.fa')
    barcode_deflines = list(barCodeDict.keys())

    info = []
    paired_fa = []
    for name in orderlist:
        print('working on %s' % name)
        fasta = '>%s(assembled)\n%s\n' % (name, seqDict[name])
        paired_fa.append(fasta)
        br_result = do_barcode_blast(fasta)
        #selfBarcode, isInBlast = find_barcode_of_itself(name.split('_')[0], br_result)
        selfBarcode = None
        if selfBarcode == None:
            #didn't find the barcode in the dataset at all
            difference, subjct = compare_top_hit(br_result)
            bfa = '>%s\n%s\n' % (rename_pair(subjct), barCodeDict[subjct])
            paired_fa.append(bfa)
            keys = list(difference.keys())
            keys.sort()
            for pos in keys:
                info.append('%s\t%s\t%s\t%s\n' % (name, subjct, pos, difference[pos]))
        else:
            #found the barcode of itself
            subjct, seq = selfBarcode
            bfa = '>%s\n%s\n' % (subjct, seq)
            difference = compare_itself(fasta, bfa)
            paired_fa.append(bfa)
            keys = list(difference.keys())
            keys.sort()
            for pos in keys:
                info.append('%s\t%s\t%s\t%s\n' % (name, subjct, pos, difference[pos]))
            if isInBlast:
                #itself is found by blast
                pass
            else:
                #this is very special
                #the barcode is not in the blast result
                difference, subjct = compare_top_hit(br_result)
                bfa = '>%s(addBack_closest_barcode)\n%s\n' % (subjct, barCodeDict[subjct])
                paired_fa.append(bfa)
                keys = list(difference.keys())
                keys.sort()
                for pos in keys:
                    info.append('%s\t%s\t%s\t%s\n' % (name, subjct, pos, difference[pos]))

        info.append('#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#\n')
    dn = cmn.lastName(fn) + '.report'
    cmn.write_file(''.join(info), dn)

    dn = cmn.lastName(fn) + '_paired.fa'
    cmn.write_file(''.join(paired_fa), dn)
