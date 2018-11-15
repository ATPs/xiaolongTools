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
def parse_protein_query(fn):
    alist = []
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue

        if line[0] == '>':
            Id = line[1:].split()[0]
            alist.append(Id)

    return alist

def parse_dna_query(fn):
    alist = []
    rset = set([])
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue

        if line[0] == '>':
            items = line[1:].strip().split()
            Id = items[0]
            direction = items[-1][:-1]
            if direction == '-':
                rset.add(Id)

            alist.append(Id)

    return alist, rset

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_blast_output(fn, rset):
    adict = {}
    for line in cmn.file2lines(fn):
        #qseqid sseqid pident evalue qlen qstart qend slen sstart send
        items = line.strip().split()
        if len(items) == 0:
            continue
        qseqid, sseqid = items[:2]
        if qseqid in adict:
            continue
        evalue = float(items[3])
        if evalue > 0.001:
            continue

        ident = float(items[2])
        if ident < 50:
            continue

        isForward = True
        sstart, send = list(map(int, items[8:10]))
        if sstart > send:
            isForward = False

        if qseqid in rset:
            isForward = (not isForward)

        adict[qseqid] = (sseqid, isForward)
    return adict

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


reverse_dict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        }

def reverse_strand(read):
    new = []
    for char in read[::-1]:
        try:
            a = reverse_dict[char]
        except KeyError:
            a = 'N'
        new.append(a)

    return ''.join(new)



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fprot, fdna = sys.argv[1:]
    except:
        print("Usage: *.py scaf, prot.out dna.out", file=sys.stderr)
        sys.exit()

    seqDict = read_fa(fn)

    ordered_proteins = parse_protein_query('/work/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/mito_scripts/ordered_Lerema_NCBI_proteins.fa')
    ordered_proteins.reverse()
    proteins = set(ordered_proteins)

    ordered_genes, reverse_genes = parse_dna_query('/work/archive/biophysics/Nick_lab/wli/project/sequencing/scripts/mito_scripts/Lerema_ordered_dna_modified.fa')
    print('reverse_genes', reverse_genes)

    prot_out = parse_blast_output(fprot, reverse_genes)
    dna_out = parse_blast_output(fdna, reverse_genes)

    missing = []
    taken = []
    for gene in ordered_genes:
        if gene in proteins:
            try:
                scaf, direction = prot_out[gene]
            except KeyError:
                missing.append(gene)
                continue
        else:
            try:
                scaf, direction = dna_out[gene]
            except KeyError:
                missing.append(gene)
                continue

        if len(taken) == 0:
            taken.append((scaf, direction))

        else:

            lastScaf, lastDirection = taken[-1]
            if scaf != lastScaf:
                taken.append((scaf, direction))
            else:
                #checking for the same scaffold
                if direction != lastDirection:
                    print('Error! direction not right for %s %s %s' % (gene, scaf, direction))


    print('scafs', taken)
    print('missing', missing)

    new = []
    a, b = 0, 0
    for scaf, isForward in taken:
        seq = seqDict[scaf]
        cov = int(scaf.split('cov')[-1].split('_')[0])
        a += cov * len(seq)
        b += len(seq)

        if not isForward:
            seq = reverse_strand(seq)
        new.append(seq)

    fasta = '>denovo_cov%s\n%s\n' % (int(a/b), 'NNNNNNNNNN'.join(new))
    dn = 'denovo4gapClose.fa'
    cmn.write_file(fasta, dn)
