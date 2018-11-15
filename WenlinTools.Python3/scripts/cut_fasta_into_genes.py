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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    orderList = []
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
        orderList.append(defline)
    return adict, orderList


def read_gene_range(fn):
    rangeDict = {}
    with open(fn) as fp:
        for line in fp:
            exon, i, j = line.strip().split()
            prot = exon.split('.e')[0]
            i, j = int(i), int(j)
            try:
                rangeDict[prot][0] = min(i, rangeDict[prot][0])
                rangeDict[prot][1] = max(j, rangeDict[prot][1])
            except KeyError:
                rangeDict[prot] = [i, j]

    return rangeDict


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, frange = sys.argv[1:]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()


    geneRange = read_gene_range(frange)

    seqDict, order_list = read_fa(fn)


    stat = []
    outdir = '%s_gene_fasta' % cmn.lastName(fn)
    cmn.mkdir(outdir)
    for gene in geneRange:
        i, j = geneRange[gene]
        print(gene, i, j)
        stat.append('%s\t%s\n' % (gene, j - i))

        dn = '%s/%s.fa' % (outdir, gene)
        with open(dn, 'w') as dp:
            for name in order_list:
                seq = seqDict[name][i:j]
                if seq.strip('-').strip('N') == '':
                    continue
                fasta = '>%s\n%s\n' % (name, seq)
                dp.write(fasta)

        if cmn.filesize(dn) == 0:
            print('fileSize0', dn)


    dn = cmn.lastName(fn) + '_takenRange.info'
    cmn.write_file(''.join(stat), dn)




