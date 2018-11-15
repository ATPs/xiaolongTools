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
import barcode_processing as bp


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

def read_baits(fn):
    adict = {}
    toAdd = {}
    hasPrimer = True
    new = []
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue
        sp, name, seq = line.split()
        print(len(seq))
        if len(seq) != 698:
            hasPrimer = False
            if len(seq) == 658:
                #fixable
                seq = add_primer(seq)
            else:
                print('Error! didn\'t recognize the length of the bait %s %s' % (sp, name))
                sys.exit()
        newline = '%s\t%s\t%s\n' % (sp, name, seq)
        new.append(newline)
        key = '%s_%s' % (sp, name)
        adict[key] = seq
        toAdd[name] = seq

    if not hasPrimer:
        print('revise the input baits to add primer...')
        cmn.write_file(''.join(new), fn)

    return adict, toAdd


def add_primer(seq):
    return 'ACTAATCATAAAGATATTGG%sTGATTTTTTGGTCATCCAGA' % seq.strip()

def parse_ref(seqDict):
    cmn.mkdir('baits')

    newDict = {}
    for i, name in enumerate(seqDict):
        seq = seqDict[name]
        fnlabel = 'bait%s' % i
        dn = 'baits/%s.fa' % fnlabel
        name = name.replace('*', '').replace('"', "'")
        fasta = '>%s\n%s\n' % (name, seq)
        cmn.write_file(fasta, dn)
        cmd = 'module add bwa; bwa index %s -p %s' % (dn, fnlabel)
        cmn.run(cmd)
        newDict[name] = dn
    return newDict

def group_fq(fqlist):
    adict = {}
    for fn in fqlist:
        fnlabel = cmn.lastName(fn)
        sp = fnlabel.split('_')[0]
        if sp not in adict:
            adict[sp] = [None, None, None]
        if '_R1' in fnlabel:
            adict[sp][0] = fn
        elif '_R2' in fnlabel:
            adict[sp][1] = fn
        elif '_singleton' in fnlabel:
            adict[sp][2] = fn
        else:
            print('Error! Can not recognize file %s' % fn)
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def old_log_newBaits_ifPossible(seqs):
    fall = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    takenNames = set([each.strip()[1:] for each in cmn.cmd2lines('grep ">" %s' % fall)])

    fnew = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/addedBaits_fromPipeline.fa'
    seqDict = read_fa(fnew)

    for name in seqs:
        if name not in takenNames and (name not in seqDict):
            seqDict[name] = seqs[name][20:678]


    with open(fnew, 'w') as dp:
        for name in seqDict:
            if name not in takenNames and (name not in seqDict):
                print('saving %s into database...' % name)
            fasta = '>%s\n%s\n' % (name, seqDict[name])
            dp.write(fasta)

    fverify = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/all_barcodes_4verify.fa'
    dict2 = read_fa(fverify)
    seqDict.update(dict2)
    with open(fverify, 'w') as dp:
        for name in seqDict:
            fasta = '>%s\n%s\n' % (name.replace('(assembled)', '').strip('.'), seqDict[name].replace('-', 'N'))
            dp.write(fasta)
    cmd = 'module add blast;cd /project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes; makeblastdb -in=all_barcodes_4verify.fa -dbtype=nucl; chmod a+w all_barcodes_4verify.*'
    cmn.run(cmd)


def log_newBaits_ifPossible(seqs):
    conn = bp.lock_database()
    for name in seqs:
        print('name', name)
        seq = seqs[name]
        if len(seq) != 658:
            seq = seq[20:678]
        items = name.strip().split('_')
        try:
            sampleID, genus, sp = items[:3]
        except:
            #if the format is not right, then it is those in the database already
            print('skip logging %s due to the format problem' % name)
            continue

        tmp, CSid = bp.find_IDs(name)
        current_list = bp.select_barcodes('(genus = \'%s\' AND species = \'%s\')' % (genus, sp), conn)
        #(sampleID text, CSid text, isPCR INTEGER, genus text, species text, full_name text PRIMARY KEY, sequence text)""")
        isInDB = False
        for each in current_list:
            dbSeq = str(each[-1])
            if seq == dbSeq:
                isInDB = True
                break
        if not isInDB:
            newName = '%s_%s|bait_%s' % (genus, sp, sampleID)
            bp.add_barcode(sampleID, 1, CSid, newName, seq, conn)
        else:
            print('%s is from the database, skip' % name)

    bp.close_database(conn)


if __name__=='__main__':
    #options=parse_options()
    try:
        fref, fqlist =sys.argv[1:3]
    except:
        print("Usage: *.py sampleInfo.baits fqlist", file=sys.stderr)
        sys.exit()

    #add primer if not added
    ref_seqs, toAddDict = read_baits(fref)
    #log the baits into the dataset
    log_newBaits_ifPossible(ref_seqs)

    #index ref here
    frefs = parse_ref(ref_seqs)

    fqlist = cmn.getid(fqlist)
    fq_groups = group_fq(fqlist)

    N = cmn.cpu_check()

    bwa_cmds = ['module add bwa']
    for reflabel in frefs:
        fref = frefs[reflabel]
        fnlabel = cmn.lastName(fref).replace('.fa', '')
        for sp in fq_groups:
            R1, R2, single = fq_groups[sp]
            cmd = 'bwa mem -t %s -B 2 -M %s %s %s | grep "%s" > %s_paired_%s_mapped.sam ' % (N, fnlabel, R1, R2, reflabel, sp, fnlabel)
            bwa_cmds.append(cmd)
            cmd = 'bwa mem -t %s -B 2 -M %s %s | grep "%s" > %s_single_%s_mapped.sam ' % (N, fnlabel, single, reflabel, sp, fnlabel)
            bwa_cmds.append(cmd)
            bwa_cmds.append('\nwait\n')

    dn = 'bwa.cmds'
    cmn.write_lines(bwa_cmds, dn)





