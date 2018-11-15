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
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_gff(fn):
    #000000F GenomeHubs      gene    187555  198145  .       +       .       ID=JC_0000001
    coding = {}
    lenDict = {}
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        if len(items) < 5:
            continue

        if items[2] == 'gene':
            i, j = list(map(int, items[3:5]))
            scaf = items[0]
            prot = items[-1].split('=')[-1]
            length = j - i
            name = '%s_%s' % (scaf, prot)
            lenDict[name] = length
            for index in range(i, j + 1):
                try:
                    coding[scaf][index] = name
                except KeyError:
                    coding[scaf] = {index:name}

    return coding, lenDict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fgff, outlabel = sys.argv[1:]
    except:
        print("Usage: *.py bam gff outlabel", file=sys.stderr)
        sys.exit()

    #{'scaf': {indexI: geneName}}
    #they are grouped by gene
    print('loading gff info...')
    gff_dict, lenDict = read_gff(fgff)

    print('reading bam file...')
    samfile = pysam.AlignmentFile(fn, 'rb')

    count_dict = {}
    for read in samfile:
        if read.is_unmapped or read.is_secondary:
            continue

        scaf = read.reference_name
        try:
            subdict = gff_dict[scaf]
        except KeyError:
            continue

        for i,j in read.aligned_pairs:
            if i == None or j == None:
                continue
            index = j + 1
            try:
                gene = subdict[index]
            except KeyError:
                continue

            if gene not in count_dict:
                count_dict[gene] = {}

            try:
                count_dict[gene][index] += 1
            except KeyError:
                count_dict[gene][index] = 1
    samfile.close()


    dn = '%s_geneCov.txt' % outlabel
    with open(dn, 'w') as dp:
        for gene in count_dict:
            Ndict = count_dict[gene]
            Nlist = list(Ndict.values())
            Nlist.sort()
            N = sum(Nlist)
            median = Nlist[len(Nlist)/2]
            length = lenDict[gene]
            cov = float(N) / length
            line = '%s\t%s\t%s\t%s\n' % (gene, length, median, cov)
            dp.write(line)



