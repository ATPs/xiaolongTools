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
    return adict, len(seq)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':

    #fn = 'all_genomes_noGap.fa'
    #fn = 'all_genomes_charGap.fa'
    try:
        fn, parts = sys.argv[1:3]
        parts = int(parts)
    except:
        print('*.py all_genomes_charGap.fa 100')
        sys.exit()

    #adict, length = read_fa(fn)

    #random sample

    #parts = 100

    fnlabel = cmn.lastName(fn).replace('.fa', '')

    outdir = 'partition_%s' % fnlabel
    cmn.run('rm -r %s 2> /dev/null' % outdir)
    cmn.mkdir(outdir)

    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].strip()
            else:
                seq = line.strip()
                length = len(seq)
                each = length / parts
                for i in range(0, length, each):
                    j = i + each
                    dn = '%s/%s_%s.fa' % (outdir, fnlabel, i)
                    subseq = seq[i:j]
                    fasta = '>%s\n%s\n' % (name, subseq)
                    cmn.append_file(fasta, dn)

