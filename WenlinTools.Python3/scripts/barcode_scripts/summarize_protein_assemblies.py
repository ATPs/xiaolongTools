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


reverse_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


def reverse_strand(read):
    new = []
    for char in read[::-1]:
        try:
            a = reverse_dict[char]
        except KeyError:
            a = 'N'
        new.append(a)

    return ''.join(new)

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def barcodeRegionMapping(qstart, qend, sstart, send, barcodeI, barcodeJ):
    #firstly, get the mapping of barcode regions
    step = 3
    if sstart > send:
        step = -3

    diff = barcodeI - qstart
    sstart += step * diff

    diff = qend - barcodeJ
    send -= step * diff
    return (sstart, send)


def make_seq(oSeq, start, end):
    seq = []
    if True:
        if start < 1:
            while(start < 1):
                seq.append('N')
                start += 1
        seq.append(oSeq[start:end])

        if end > len(oSeq):
            while(end > len(oSeq)):
                seq.append('N')
                end -= 1
    return ''.join(seq)


def make_barcode(oSeq, xxx_todo_changeme):
    (start, end) = xxx_todo_changeme
    seq = []
    if start < end:
        #forward
        seq = make_seq(oSeq, start, end)

    else:
        #reverse
        seq = make_seq(oSeq, end, start)
        seq = reverse_strand(seq)

    seq = seq[1:]
    return seq




if __name__=='__main__':
    wdir = sys.argv[1]
    os.chdir(wdir)

    fns = cmn.cmd2lines('ls */*cox1_contig.fa')

    dp = open('all_prot_contigs.fa', 'w')
    for fn in fns:
        sample = cmn.lastName(fn).split('.cox1')[0]
        seqDict = read_fa(fn)
        if len(seqDict) == 1:
            fasta = '>%s\n%s\n' % (sample, list(seqDict.values())[0])
            dp.write(fasta)

        else:
            for name in seqDict:
                fasta = '>%s_%s\n%s\n' % (sample, name, seqDict[name])
                dp.write(fasta)

    dp.close()

    fass = 'all_prot_contigs.fa'
    cmd = 'makeblastdb -in=%s -dbtype=nucl; tblastn -query /project/biophysics/Nick_lab/wli/sequencing/scripts/barcode_scripts/cox1.fa -db %s ' % (fass, fass)
    cmd += '-out cox1_to_ass.out -outfmt \'6 qseqid sseqid qlen qstart qend slen sstart send\''
    cmn.run(cmd)

    barcodeI, barcodeJ = 13, 233
    mapDict = {}
    for line in cmn.file2lines('cox1_to_ass.out'):
        items = line.strip().split()
        qseqid, sseqid, qlen, qstart, qend, slen, sstart, send = items
        qstart = int(qstart)
        qend = int(qend)
        if qstart > barcodeJ or qend < barcodeI:
            continue
        else:
            mapDict[sseqid] = barcodeRegionMapping(int(qstart), int(qend), int(sstart), int(send), barcodeI, barcodeJ)


    seqDict = read_fa('all_prot_contigs.fa')

    takenID = set([])
    dp = open('all_protBarcodes.fa', 'w')
    dp_good = open('all_protBarcodes_complete.fa', 'w')
    appearing = {}
    for name in seqDict:
        try:
            mapping = mapDict[name]
        except KeyError:
            continue

        print(name, mapping)
        seq = make_barcode(seqDict[name], mapping)
        if mapping[0] > mapping[1]:
            print('reversed sequence %s' % name)
            tSeq = seq[1:-2]
        else:
            print('forward sequence %s' % name)
            tSeq = seq[:-3]
        fasta = '>%s\n%s\n' % (name, tSeq)
        dp.write(fasta)
        if 'N' not in seq and len(seq) == 661:
        #if 'N' not in seq:
            takenID.add(name.split('_')[0])
            dp_good.write(fasta)

    dp_good.close()
    dp.close()


