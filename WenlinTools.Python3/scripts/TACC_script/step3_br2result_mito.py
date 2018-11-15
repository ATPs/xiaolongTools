#!/usr/bin/python
import os, sys, numpy, string

python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
from collections import Counter

reverse_dict = {
    'A': 'T',
    'T': 'A',
    'G': 'C',
    'C': 'G'
}

def reverse_strand(seq):
    global reverse_dict
    new = []
    for char in seq[::-1]:
        try:
            rChar = reverse_dict[char]
        except KeyError:
            rChar = 'N'
        new.append(rChar)

    new = ''.join(new)
    return new

def diff_letters(a,b):
    return sum (a[i] != b[i] for i in range(len(a)))

def most_common(list):
	counts = {}
	for item in list:
		try:
			counts[item] += 1
		except KeyError:
			counts[item] = 1
	for item in list(counts.keys()):
		if counts[item] >= 0.8 * len(list) and len(list) >= 4:
			return item 
	return "-"

def most_common_biallelic(list):
    counts = {}
    for item in list:
        try:
            counts[item] += 1
        except KeyError:
            counts[item] = 1
    items = []
    for item in list(counts.keys()):
        if counts[item] >= 0.8 * len(list):
            items.append(item)
            items.append(item)
        elif counts[item] >= 0.4 * len(list):
            items.append(item)
    

    if len(items) == 2:
        if items[0] != items[1]:
            print('warning: heterozyocity position found!', file=sys.stderr)
        return items
    else:
        return ['-','-']


def tell_dominate_char(chars):
    count_dict = Counter(chars)
    print(count_dict)
    mainChar = max(list(count_dict.keys()), key=lambda x: count_dict[x])
    N = count_dict[mainChar]
    percent = float(N) / len(chars)

    if percent > 0.8 and N > 4:
        return mainChar
    
    return 'N'


try:
    fblast, fread = sys.argv[1:3]
except:
    print('*.py blastExon.dict.pkl readSp.dict.pkl', file=sys.stderr)
    sys.exit()


blast_dict = cmn.pickle_read(fblast)

#reads[sp][name] = seq
read_dict = cmn.pickle_read(fread)
splist = set(read_dict.keys())


#read in ranges for exons
#{'COX1': readID: {2:(1,2,3)}}
#TODO: fill in missing

stack_dict = {}
exon_lengths = {}
for exon in blast_dict:
    info = blast_dict[exon]
    stack_dict[exon] = {}

    #1. get the ID list
    exon_reads = {line.split()[2]: line for line in info}

    for sp in read_dict:
        stack_dict[exon][sp] = {}
        sp_reads = read_dict[sp]

        readIDs = set(sp_reads.keys()) & set(exon_reads.keys())
        
        if len(readIDs) == 0:
            print('missing regions for %s in %s' % (exon, sp), file=sys.stderr)
            continue

        for readID in readIDs:
            seq = sp_reads[readID]
            #sp qseqid sseqid pident evalue qlen qstart qend slen sstart send qseq sseq
            blastHit = exon_reads[readID]
            items = blastHit.strip().split()
            sid = items[2]
            qlen, qstart, qend, slen, sstart, send = list(map(int, items[5:11]))
            exon_lengths[exon] = qlen
        
        
            #get the corresponding region of hit
            if sstart < send:#forward strand
                hsp = seq[sstart-1: send]
                shsp_length = send - sstart + 1
            else: # backward strand
                hsp = seq[send - 1: sstart]
                hsp = reverse_strand(hsp)
                shsp_length = sstart - send + 1
            
            qhsp_length = qend - qstart + 1
            #print qstart, qend, sstart, send, qhsp_length, shsp_length
            if qhsp_length *3 != shsp_length:
                print('detected gap for %s, skip' % sid, file=sys.stderr)
                continue
            
            qseq, sseq = items[-2:]
            if ('-' in qseq) or ('-' in sseq):
                print('detected gap for %s, skip' % sid, file=sys.stderr)
                continue
            
            for i in range(qhsp_length):
                qI = i + qstart
                sI = 3 * i
                sJ = sI + 3
                codon = hsp[sI: sJ]
                try:
                    stack_dict[exon][sp][qI].append(codon)
                except KeyError:
                    stack_dict[exon][sp][qI] = [codon]



            
#tell and output
dn = 'assembled.contigs'
new = []
cov_info = []
for exon in stack_dict:
    length = exon_lengths[exon]
    spDict = stack_dict[exon]
    for sp in spDict:
        stacks = spDict[sp]

        newSeq = []

        cov = 0

        for each in range(length):
            p = each + 1
            print(p, exon, sp)
            try:
                alist = stacks[p]
                N = len(alist)
                cov += N
            except:
                newSeq.append('---')
                continue

            for i in range(3):
                chars = [each[i] for each in alist]
                main_char = tell_dominate_char(chars)
                newSeq.append(main_char)
        cov = float(cov) / length
        cov_info.append('%s\t%s\t%s\n' % (exon, sp, cov))

        line = '%s\t%s\t%s\n' % (exon, sp, ''.join(newSeq))
        new.append(line)

cmn.write_file(''.join(new), dn)
cmn.write_file(''.join(cov_info), 'cov_info.txt')


