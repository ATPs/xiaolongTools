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
    for line in cmn.file2lines(fn):
        sp, name, seq = line.split()
        key = '%s_%s' % (sp, name)
        adict[key] = seq
    return adict


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

def add_primer(seq):
    return 'ACTAATCATAAAGATATTGG%sTGATTTTTTGGTCATCCAGA' % seq.strip()


def add_indel(seq, indel_dict):
    seq = list(seq)
    for i, j in indel_dict:
        indel = 'N' * len(indel_dict[(i,j)])
        seq[i] += indel
    return ''.join(seq)

def insert_in_ref_info(info, indel_dict):
    fasta = info.split('>')[1:]
    new = []
    for each in fasta:
        name, seq = each.strip().split('\n')
        seq = add_indel(seq, indel_dict)
        fasta = '>%s\n%s\n' % (name, seq)
        new.append(fasta)
    return ''.join(new)


def insert_in_lines(lines, indel_dict):
    new = []
    for line in lines:
        sp, defline, seq = line.strip().split()
        seq = add_indel(seq, indel_dict)
        new.append('%s\t%s\t%s' % (sp, defline, seq))
    return new


def read_indel_info(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        key, a, b, char = line.strip().split('\t')
        i, j = list(map(int, key[1:-1].split(', ')))
        adict[(i,j)] = char
    return adict


def add_in_baits(fref):
    fbait = 'sampleInfo.baits'
    ref_info = cmn.txt_read(fref)
    if cmn.filexist('bait_insertion'):
        indel_dict = read_indel_info('bait_insertion')
    else:
        indel_dict = {}

    #baits added by customBaits
    #fadd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/added_from_customBaits.baitInfo'
    #if cmn.filexist(fadd):
    #    add_lines = cmn.file2lines(fadd)
    #else:
    #    add_lines = []
    add_lines = []

    if len(indel_dict) != 0:
        ref_info = insert_in_ref_info(ref_info, indel_dict)
        add_lines = insert_in_lines(add_lines, indel_dict)



    refIDs = [line[1:] for line in ref_info.split('\n')
            if line.strip() != '' and line[0] == '>']
    addedIDs = [line.split()[1] for line in add_lines]

    #new = []
    #when check a new line, need to check both the fref and the fadd
    #if the one is not in fadd, add it to fadd
    for line in cmn.file2lines(fbait):
        sp, defline, seq = line.strip().split()
        if all([defline.upper() not in refID.upper() for refID in refIDs]):
            #not in ref
            if all([defline.upper() not in addedID.upper() for addedID in addedIDs]):
                if len(seq) == 698:
                    add_lines.append(line)
                else:
                    if len(seq) != 658:
                        print('Error! length of bait barcode is wrong for %s %s' % (sp, defline))
                        sys.exit()
                    else:
                        seq = add_primer(seq)
                        add_lines.append('%s\t%s\t%s' % (sp, defline, seq))

    #now get a new fadd, need to format it into fasta
    add_fasta = []
    for line in add_lines:
        sp, defline, seq = line.strip().split()
        fasta = '>%s\n%s\n' % (defline, seq)
        add_fasta.append(fasta)

    ref_info += '\n'
    ref_info += ''.join(add_fasta)
    dn = 'species_barcodes_4mapping_withAddon.fa'
    cmn.write_file(ref_info, dn)

    #index it
    cmd = 'module add bwa; bwa index %s' % dn
    cmn.run(cmd)

    #record the new fadd
    #cmn.write_lines(add_lines, fadd)

    return dn



if __name__=='__main__':
    #options=parse_options()
    try:
        fqlist =sys.argv[1]
    except:
        print("Usage: *.py fqlist", file=sys.stderr)
        sys.exit()

    #ref_seqs = read_baits(fref)
    #index ref here
    #frefs = parse_ref(ref_seqs)

    #need to read in current baits to put into the libs
    fref = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    #assume the bait file name is sampleInfo.baits, add it in the libs
    fref = add_in_baits(fref)
    #cmn.run('ls %s*| awk \'{printf "ln -s %%s\\n", $1}\'|bash' % fref)
    frefs = {'barcode': cmn.lastName(fref)}

    fqlist = cmn.getid(fqlist)
    fq_groups = group_fq(fqlist)

    N = cmn.cpu_check()
    bwa_cmds = ['module add bwa']
    for reflabel in frefs:
        fref = frefs[reflabel]
        fnlabel = fref
        for sp in fq_groups:
            R1, R2, single = fq_groups[sp]
            cmd = 'bwa mem -B 2 -t %s -M %s %s %s > %s_paired_%s.sam ' % (N, fnlabel, R1, R2, sp, fnlabel)
            bwa_cmds.append(cmd)
            cmd = 'bwa mem -B 2 -t %s -M %s %s > %s_single_%s.sam ' % (N, fnlabel, single, sp, fnlabel)
            bwa_cmds.append(cmd)
            #bwa_cmds.append('\nwait\n')

    dn = 'bwa2.cmds'
    cmn.write_lines(bwa_cmds, dn)





