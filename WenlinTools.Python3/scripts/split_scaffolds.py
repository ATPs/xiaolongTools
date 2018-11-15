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
    info = cmn.txt_read(fa).decode("utf-8", "ignore")
    fastas = info.split('>')[1:]
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

    #fn = 'all_genomes_noGap.fa'
    #fn = 'all_genomes_charGap.fa'
    try:
        fn = sys.argv[1]
    except:
        print('*.py all_genomes_charGap.fa ')
        sys.exit()

    adict = read_fa(fn)


    fnlabel = cmn.lastName(fn).replace('.fa', '')
    outdir = 'splitS_%s' % fnlabel
    cmn.mkdir(outdir)
    for i, key in enumerate(adict):
        seq = adict[key]
        fasta = '>%s\n%s\n' % (key, seq)
        dn = '%s/%s_%s.fa' % (outdir, fnlabel, i)
        cmn.write_file(fasta, dn)

