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
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    fn1 = 'sum_denovo.fa'
    fn2 = 'sum_barcodes.fa'
    #fn3 = 'compare.check'

    if 'Error' in cmn.txt_read('compare.check'):
        print('Error in running barcode pipeline! please fix lines with "Error" in "compare.check" file!')
        sys.exit()

    replaceIDs = set(cmn.cmd2lines('grep takenD compare.check|grep -v same|cut -f 1'))

    seqDict1 = read_fa(fn1)
    seqDict2 = read_fa(fn2)

    dn = 'sum_hybrid.fa'
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
                            print('please manually fix %s (length:%s)' % (name, len(seq)))

            else:
                seq = seqDict2[name]

            fasta = '>%s\n%s\n' % (name, seq)
            fp.write(fasta)

