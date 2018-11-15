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
rdict = {
        'T': 'A',
        'A': 'T',
        'G': 'C',
        'C': 'G',
        }

def reverse_sequence(seq):
    new = []
    for char in seq[::-1]:
        try:
            new.append(rdict[char])
        except KeyError:
            new.append('N')
    return ''.join(new)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fn):
    fasta = cmn.txt_read(fn).strip().split('>')[1:]

    seqDict = {}
    orderlist = []
    for each in fasta:
        lines = each.split('\n')
        defline = lines[0].split()[0]
        seq = ''.join(lines[1:])
        seqDict[defline] = seq
        orderlist.append(defline)
    return seqDict, orderlist

def read_mito(fn):
    seq = [line.strip() for line in cmn.file2lines(fn)
            if line[0] != '>']
    return ''.join(seq)


def parse_blastDict(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        #qseqid sseqid evalue pident qlen qstart qend slen sstart send
        items = line.strip().split()
        gene = items[0]
        qlen, qstart, qend, slen, sstart, send = list(map(int, items[4:10]))
        i, j, isReverse = find_start_and_end(qlen, qstart, qend, slen, sstart, send)
        adict[gene] = (i, j, isReverse)
    return adict


def find_start_and_end(qlen, qstart, qend, slen, sstart, send):
    if sstart < send:
        #forward strand
        if qstart != 1:
            sstart -= ( qstart - 1 ) * 3
        if qend != qlen:
            send += ( qlen - qend ) * 3
        isReverse = False
        data = (sstart, send, isReverse)
    else:
        if qstart != 1:
            sstart += ( qstart - 1 ) * 3
        if qend != qlen:
            send -= ( qlen - qend ) * 3
        isReverse = True
        data = (send, sstart, isReverse)
    return data

if __name__=='__main__':
    #options=parse_options()
    try:
        fprot, fmito, fblast = sys.argv[1:]
    except:
        print("Usage: *.py query_protein.fa mito.fa blast.out", file=sys.stderr)
        sys.exit()

    seqDict, order_list = read_fa(fprot)

    order_list = [line.split()[0] for line in cmn.file2lines('dnaRange.txt')]

    baseSeq = read_mito(fmito)

    blastDict = parse_blastDict(fblast)

    seqlist = []
    for gene in order_list:
        try:
            i, j, isReverse = blastDict[gene]
            seg = baseSeq[i-1:j]
            if isReverse:
                seg = reverse_sequence(seg)
            print(gene, j-i+1, i, j, isReverse)
        except KeyError:
            seg = 'N' * len(seqDict[gene] * 3)
            print(gene, len(seg), 'notFound')

        seqlist.append(seg)

    dn = 'algined.fa'
    fasta = '>aligned\n%s\n' % (''.join(seqlist))
    cmn.write_file(fasta, dn)

