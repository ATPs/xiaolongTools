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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:]).strip()
        adict[defline] = seq
    return adict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def parse_IDmapping_and_newDict(fn):
    print('processing fasta names...')
    cmd = 'source /home2/wli/.bash_profile; rename_fa_fullname.py %s > tmp.namelog' % fn
    cmn.run(cmd)

    seqDict = read_fa(fn + '.renamed')
    IDmapping = names2IDs(list(seqDict.keys()))
    return IDmapping, seqDict

def names2IDs(alist):
    adict = {}
    for name in alist:
        ID = name.split('_')[0].split(']')[-1]
        adict[ID] = name
    return adict


def load_verified_barcodes():
    fgood = '/archive/biophysics/Nick_lab/wli/archive/barcodes/auto_tables/verified_barcodes.fa'

    seqDict = read_fa(fgood)
    IDmapping = names2IDs(list(seqDict.keys()))
    return IDmapping, seqDict



if __name__=='__main__':
    fn1 = 'sum_denovo.fa'
    fn2 = 'sum_barcodes.fa'
    #fn3 = 'compare.check'

    if 'Error' in cmn.txt_read('compare.check'):
        print('##########################################################################')
        print('Error in running barcode pipeline! please fix lines with "Error" in "compare.check" file!')
        print('##########################################################################')
        #sys.exit()

    replaceIDs = set(cmn.cmd2lines('grep takenD compare.check|grep -v same|cut -f 1'))

    seqDict1 = read_fa(fn1)
    seqDict2 = read_fa(fn2)

    dn = 'sum_hybrid.fa'
    #newDict = {}
    with open(dn, 'w') as fp:
        for name in seqDict2:
            if name in replaceIDs:
                seq = seqDict1[name]
                name = name + '_denovo'
                if len(seq) != 658:
                    diffN = len(seq) - 658
                    if seq.count('N') == diffN:
                        seq = seq.replace('N', '')
                    else:
                        if seq[:diffN].count('N') == diffN:
                            seq = seq[diffN:]
                        else:
                            print('Error: please manually fix %s (length:%s)' % (name, len(seq)))

            else:
                seq = seqDict2[name]

            #newDict[name] = seq
            if len(seq) != 658:
                seq = seq[:658]
                addN = 658 - len(seq)
                if addN > 0:
                    seq += '-' * addN
            fasta = '>%s\n%s\n' % (name, seq)
            fp.write(fasta)

    ####backup sum_hybrid.fa

    #1. get IDs of the sequences
    IDmapping, newDict = parse_IDmapping_and_newDict(dn)

    #2. load done sequences
    doneIDdict, doneDict = load_verified_barcodes()

    #3. construct the sequences that put into intermediate
    backup_fasta = []
    for ID in IDmapping:
        label = ''
        if ID in doneIDdict:
            #label = ''
            name = doneIDdict[ID]
            seq = doneDict[name]
        else:
            #label = '[pp]'
            name = IDmapping[ID]
            seq = newDict[name]

        fasta = '>%s%s\n%s\n' % (label, name, seq)
        backup_fasta.append(fasta)

    dn = 'fasta4backup.fa'
    cmn.write_file(''.join(backup_fasta), dn)

    #make the backup
    remote_dn = '/data/www/wenlin/html/transfer/intermediate_barcodes/tmp4backup.fa'
    cmd = 'rsync %s wenlin@morpho.swmed.edu:%s' % (dn, remote_dn)
    print('transfering the file to morpho server...')
    cmn.run(cmd)

    inCmd = 'python /home/wenlin/my_programs/backup_intermediate_barcodes.py %s' % remote_dn
    cmd = 'ssh -t wenlin@morpho.swmed.edu "%s"' %  inCmd
    print('running remote command...')
    cmn.run(cmd)
