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
def get_query_sequence(seqDict, genus, sp):
    #1. anything in Eudamine file has higher priority
    fEud = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/Eudaminae-barcode-reference.txt'
    cmd = 'grep %s %s' % (sp, fEud)
    lines = cmn.cmd2lines(cmd)
    if len(lines) == 1:
        name = lines[0].split()[0]
        seq = seqDict[name]
        fasta = '>%s\n%s\n' % (name, seq)
        qlen = len(seq.replace('N', ''))
        print('pick %s for %s %s' % (name, genus, sp))
        return fasta, qlen

    #look it up in other files
    names = list(seqDict.keys())
    good_names = [name for name in names
            if genus in name ]
    if len(good_names) == 0:#sp is just 'sp'
        print('can not find barcode for genus keyword "%s"' % genus)
        good_names = names

    if len(good_names) > 1:
        #try to refine it
        tmp = [name for name in good_names
                if sp in name]
        if len(tmp) != 0:
            good_names = tmp

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
    seq = seqDict[name]
    fasta = '>%s\n%s\n' % (name, seq)
    qlen = len(seq.replace('N', ''))
    print('pick %s for %s %s' % (name, genus, sp))
    return fasta, qlen

def do_barcode_blast(sequence):
    fdb = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #fdb = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_NoN_0.95.fasta'
    namelabel = sequence.split()[0][1:].split()[0].split('|')[0].replace('*','').split('[')[0].replace('"','').replace("'", '')
    fquery = '/tmp/%s.fa' % namelabel
    cmn.write_file(sequence, fquery)
    cmd = 'module add blast; blastn -query %s -db %s ' % (fquery, fdb)
    cmd += '-outfmt \'6 sseqid qlen slen length pident\''
    lines = cmn.cmd2lines(cmd)
    #cmn.run('rm %s' % fquery)
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
    Nbait = 10

    baits = []

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
                sseq = seqDict[sseqid]
                slen = len(sseq.replace('N', ''))
                if slen < qlen:
                    continue

                #if reach here, everything is good
                print('taken bait: %s %s %s' % (ident, slen, sseqid))
                #fasta = '>%s\n%s\n' % (sseqid, sseq)
                baits.append((sseqid, sseq))
        identCut -= 1
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

    info = []
    for line in cmn.file2lines(fn):
        #5077    Autochton zarex
        items = line.strip().split()
        sample, genus, sp = items[:3]
        query_sequence, qlen = get_query_sequence(seqDict, genus, sp)
        br_result = do_barcode_blast(query_sequence)
        print(br_result)
        #print '\n'.join(br_result)
        baits = pick_barcode_baits(br_result, qlen, seqDict)
        info += format_baits(sample, baits)

    dn = cmn.lastName(fn) + '.baits'
    cmn.write_file(''.join(info), dn)



