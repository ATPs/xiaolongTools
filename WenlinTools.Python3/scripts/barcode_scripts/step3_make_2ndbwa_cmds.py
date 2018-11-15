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
    fref = '/project/biophysics/Nick_lab/wli/sequencing/scripts/data/barcodes/species_barcodes_4mapping.fa'
    cmn.run('ls %s*| awk \'{printf "ln -s %%s\\n", $1}\'|bash' % fref)
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
            cmd = 'bwa mem -t %s -M %s %s %s > %s_paired_%s.sam ' % (N, fnlabel, R1, R2, sp, fnlabel)
            bwa_cmds.append(cmd)
            cmd = 'bwa mem -t %s -M %s %s > %s_single_%s.sam ' % (N, fnlabel, single, sp, fnlabel)
            bwa_cmds.append(cmd)
            #bwa_cmds.append('\nwait\n')

    dn = 'bwa2.cmds'
    cmn.write_lines(bwa_cmds, dn)





