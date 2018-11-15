#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

def parse_name(name):
    name = name.replace('Pterourus', 'Papilio')
    #name = name.replace('Arzecla', 'Arumecla')
    #name = name.replace('Badecla', 'Arumecla')
    name = name.replace('(assembled)', '')
    name = name.replace(' ', '_').replace('__', '_')
    return name

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = parse_name(lines[0])
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_query_sequence(seqDict, genus, sp):
    #1. anything in Eudamine file has higher priority
    #fEud = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/Eudaminae-barcode-reference.txt'
    #cmd = 'grep %s %s' % (sp, fEud)
    #lines = cmn.cmd2lines(cmd)
    #if len(lines) == 1:
    #    name = lines[0].split()[0]
    #    seq = seqDict[name]
    #    fasta = '>%s\n%s\n' % (name, seq)
    #    qlen = len(seq.replace('N', ''))
    #    print 'pick %s for %s %s' % (name, genus, sp)
    #    return fasta, qlen

    names = list(seqDict.keys())
    #try to look up the exact match first
    expected_name = '%s_%s' % (genus, sp)
    tmp = [name for name in names
        if name.upper() == expected_name.upper()]

    if len(tmp) != 0:
        name = tmp[0]
        print('found exact match %s' % name)
        seq = seqDict[name]
        fasta = '>%s\n%s\n' % (name, seq)
        qlen = len(seq.replace('N', ''))
        return fasta, qlen


    #look it up in other files
    good_names = [name for name in names
    #        if genus.upper() in name.upper().split('_')]
            if genus.upper() == name.upper().split('_')[0]]

    useGenus = False
    if len(good_names) > 0:
        useGenus = True

    cmn.run('rm pickingLog.txt 2> /dev/null')
    if len(good_names) == 0:#sp is just 'sp'
        print('can not find barcode for genus keyword "%s"' % genus)
        good_names = names
        cmn.write_file('noGenus\n', 'pickingLog.txt')

    if len(good_names) > 1:
        #try to refine it
        tmp = [name for name in good_names
                if sp.upper() in name.upper().split('_')]
        if len(tmp) != 0:
            good_names = tmp
        else:
            cmn.append_file('noSpecies\n', 'pickingLog.txt')

    #############################################
    ####new here, auto pick sequences for those has no info
    #############################################
    if cmn.filexist('pickingLog.txt'):
        print('automatically pick bait by fastq similarity')
        fsp = 'restricted_genus.info'
        if useGenus and (not cmn.filexist(fsp)):
            cmd = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/auto_rebait.py fqlist %s' % genus
        else:
            cmd = '/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/barcode_scripts/auto_rebait.py fqlist '
        cmn.run(cmd)
        good_names = cmn.file2lines('picked_bait.txt')
        cmn.write_file('pickClosed\n', 'pickingLog.txt')


    #############################################
    #############################################
    #############################################

    #try to see if type species is there
    tmp = [name for name in good_names
            if name[0] == '*']
    if len(tmp) != 0:
        good_names = tmp
    else:
        tmp = [name for name in good_names
            if '*' in name]
        if len(tmp) != 0:
            good_names = tmp

    #then randomly pick one, get the max length ones
    name = max(good_names, key=lambda x: len(seqDict[x].replace('N', '-')))
    #name = name.replace('/', '_')
    seq = seqDict[name]
    fasta = '>%s\n%s\n' % (name, seq)
    qlen = len(seq.replace('N', ''))
    print('pick %s for %s %s' % (name, genus, sp))
    return fasta, qlen

def makeBlastDatabase(seqDict):
    dn = 'db4picking.fa'
    new = ['>%s\n%s\n' % (name, seqDict[name])
        for name in seqDict
		if seqDict[name].strip('N-X') != '']
    cmn.write_file(''.join(new), dn)
    cmd = 'module add blast; makeblastdb -dbtype=nucl -in=%s' % dn
    cmn.run(cmd)
    return dn


def do_barcode_blast(sequence, seqDict):
    #fref = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #fadd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/added_from_customBaits.baitInfo'

    fdb = makeBlastDatabase(seqDict)

    #fdb = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_NoN_0.95.fasta'
    namelabel = sequence.split()[0][1:].split()[0].split('|')[0].replace('*','').split('[')[0].replace('"','').replace("'", '')
    namelabel = namelabel.replace('/', '_')
    fquery = '/tmp/%s.fa' % namelabel
    cmn.write_file(sequence, fquery)
    cmd = 'module add blast; blastn -query %s -db %s ' % (fquery, fdb)
    cmd += '-outfmt \'6 sseqid qlen slen length pident\''
    lines = cmn.cmd2lines(cmd)
    cmn.run('rm %s' % fquery)
    return lines


def pick_barcode_baits(lines, qlen, seqDict):
    #bin the lines by length
    lenDict = {}
    for line in lines:
        slen = int(line.split()[2])
        try:
            lenDict[slen].append(line)
        except:
            lenDict[slen] = [line]

    keys = list(lenDict.keys())
    keys.sort()
    keys.reverse()

    identCut = 100 #the identity variable used to scan the result
    Nbait = 5
    lower_identCut = 95 # the lowest identity to make sure they are in a same family
    lower_identCut2 = 90 # if bait is only one, lower the cutoff

    baits = []
    backup_baits = []

    while(len(baits) < Nbait):

        for key in keys:
            lines = lenDict[key]
            for line in lines:
                if len(baits) >= Nbait:
                    break

                items = line.strip().split()
                #filtering identity
                ident = int(round(float(items[-1]), 0))
                if ident != identCut:
                    continue
                #filtering length
                sseqid = items[0]
                try:
                    sseq = seqDict[sseqid]
                except KeyError:
                    sseq = seqDict[sseqid + '.']

                #slen = len(sseq.replace('N', ''))
                #if slen < qlen:
                #    continue

                #if reach here, everything is good
                #fasta = '>%s\n%s\n' % (sseqid, sseq)
                if identCut >= lower_identCut:
                    print('taken bait: %s %s %s' % (ident, slen, sseqid))
                    baits.append((sseqid, sseq))
                if identCut >= lower_identCut2:
                    backup_baits.append((sseqid, sseq))
        identCut -= 1
        if identCut < lower_identCut2:
            break

    #if we only have one bait, take the next one with identity >= 90
    if len(baits) == 1:
        print('lower cutoff from %s to %s' % (lower_identCut, lower_identCut2))
        baits = backup_baits[:3]

    return baits



def format_baits(sp, baits):
    lines = []
    for bait in baits:
        line = '%s\t%s\t%s\n' % (sp, bait[0].replace(' ', '_'), bait[1])
        lines.append(line)
    return lines

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_baitInfo(fadd):
    adict = {}
    for line in cmn.file2lines(fadd):
        sp, defline, seq = line.split()
        defline = parse_name(defline)
        adict[defline] = seq
    return adict


def read_autoTable(fn):
    seqDict = read_fa(fn)

    new = {}
    for name in seqDict:
        seq = seqDict[name]
        name = '_'.join(name.split('_')[1:])
        new[name] = seq
    return new


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes.fasta'
    seqDict = read_fa(fall)
    #fadd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/addedBaits_fromPipeline.fa'
    #if cmn.filexist(fadd):
    #    seqDict.update(read_fa(fadd))
    #ftable = '/archive/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/verified_barcodes.fa'
    #seqDict.update(read_autoTable(ftable))

    info = []
    for line in cmn.file2lines(fn):
        #5077    Autochton zarex
        items = line.strip().replace('?', ' ').split()
        try:
            sample, genus, sp = items[:3]
        except:
            sample, genus = items[:2]
            sp = 'sp'
        sp = sp.split('-')[0].split('_')[0]
        genus = parse_name(genus)
        query_sequence, qlen = get_query_sequence(seqDict, genus, sp)
        br_result = do_barcode_blast(query_sequence, seqDict)
        print('\n'.join(br_result))
        baits = pick_barcode_baits(br_result, qlen, seqDict)
        info += format_baits(sample, baits)

    dn = cmn.lastName(fn) + '.baits'
    cmn.write_file(''.join(info), dn)



