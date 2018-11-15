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

def tell_dominate_char_2phase(chars, cov):
    count_dict = Counter(chars)
    
    #sorted from the smallest to biggest
    alist = sorted(list(count_dict.keys()), key=lambda x: count_dict[x])
    #the cutoff to claim multiple phases
    #this need to be stringent
    phaseCut = 30 #for each
    singleCut = 4
    cutoff =  max(cov * 0.20, phaseCut)

    flist = [char for char in alist
        if count_dict[char] > cutoff and (count_dict[char] > (0.2 * len(chars)))]
    
    print(count_dict)
    if len(flist) == 0:
        #coverage is low, take the most dominated one
        char1 = alist[-1]
        N = count_dict[char1]
        #if coverage is extremly low, ignore it
        if N <= singleCut:
            return ['N', 'N']

        percent = float(N) / len(chars)
        if percent > 0.8 and N > singleCut:
            return [char1, char1]
        else:
            #try to get two chase
            two_chars = alist[-2:]
            isGood = True
            for char in two_chars:
                N = count_dict[char]
                percent = float(N) / len(chars)
                print(char, N, len(chars), percent)
                if percent < 0.4 or N < phaseCut:
                    isGood = False
            if isGood:
                print(count_dict, two_chars, count_dict, cutoff, 'twochars')
                #print two_chars
                return two_chars
            else:
                return ['N', 'N']
    elif len(flist) == 1:
        char1 = flist[0]
        return [char1, char1]
    else:
        char1, char2 = flist[-2:]
        print(count_dict, char1, char2, count_dict, cutoff, 'twochars')
        return [char1, char2]

def consensus_seqs(alist):
    seq = []
    k = max([len(each) for each in alist])
    for i in range(k):
        chars = []
        for each in alist:
            try:
                char = each[i]
                if char.strip() != '':
                    chars.append(char)
            except:
                continue
        
        count_dict = Counter(chars)
        maxChar = max(list(count_dict.keys()), key=lambda x: count_dict[x])
        N = count_dict[maxChar]
        #print 'ccSee', count_dict
        if N <= 10:
            seq.append('N')
        else:
            fraction = float(N) / len(chars)
            if fraction >= 0.8:
                seq.append(maxChar)
            else:
                seq.append('N')
    return ''.join(seq)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

try:
    fblast, fread = sys.argv[1:3]
except:
    print('*.py blastExon.dict.pkl readSp.dict.pkl', file=sys.stderr)
    sys.exit()

extend = 50#extend 50bp on each side
query_length = 658

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
        end_seqs = []
        start_seqs = []
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
                #a. get the ending index corresponding to query
                endIndex = len(seq) - send + qend
                if endIndex > query_length:
                    end_seqs.append(seq[query_length-endIndex:])
                startIndex = sstart - qstart
                if startIndex > 0:
                    start_seqs.append(seq[:startIndex])

            else: # backward strand
                hsp = seq[send - 1: sstart]
                hsp = reverse_strand(hsp)
                seq = reverse_strand(seq)
                shsp_length = sstart - send + 1
                #a. get the ending index corresponding to query
                send, sstart = sstart, send
                endIndex = len(seq) - send + qend
                if endIndex > query_length:
                    end_seqs.append(seq[query_length-endIndex:])
                startIndex = sstart - qstart
                if startIndex > 0:
                    start_seqs.append(seq[:startIndex])
            
            qhsp_length = qend - qstart + 1
            #print qstart, qend, sstart, send, qhsp_length, shsp_length
            if qhsp_length != shsp_length:
                print('detected gap for %s, skip' % sid, file=sys.stderr)
                continue
            
            qseq, sseq = items[-2:]
            if ('-' in qseq) or ('-' in sseq):
                print('detected gap for %s, skip' % sid, file=sys.stderr)
                continue
            
            for i in range(qhsp_length):
                qI = i + qstart
                sI = i
                codon = hsp[sI]
                try:
                    stack_dict[exon][sp][qI].append(codon)
                except KeyError:
                    stack_dict[exon][sp][qI] = [codon]

        dn = '%s_%s_endExtend.fa' % (sp, exon)
        cmn.write_lines(end_seqs, dn)
        
        print('endseq', consensus_seqs(end_seqs))
        maxLength = max([len(seq) for seq in start_seqs])
        pattern = '{:>%s}' % maxLength
        start_seqs = [pattern.format(seq) for seq in start_seqs]
        dn = '%s_%s_startExtend.fa' % (sp, exon)
        cmn.write_lines(start_seqs, dn)
        print('startseq', consensus_seqs(start_seqs))


            
#tell and output
dn = 'phased_assemblies.contigs'
new = []
cov_dict = {}
#get coverage first
for exon in stack_dict:
    cov_dict[exon] = {}
    for sp in stack_dict[exon]:
        stacks = stack_dict[exon][sp]
        Nlist = [len(stacks[key]) for key in stacks]
        cov = float(sum(Nlist)) / len(Nlist)
        cov_dict[exon][sp] = cov

cov_info = ['%s\t%s\t%s\n' % (sp, exon, cov_dict[exon][sp])
            for exon in cov_dict for sp in cov_dict[exon]]

cmn.write_file(''.join(cov_info), 'cov_info.txt')


for exon in stack_dict:
    length = exon_lengths[exon]
    spDict = stack_dict[exon]
    for sp in spDict:
        stacks = spDict[sp]

        newSeq = [[], []]

        cov = cov_dict[exon][sp]

        for each in range(length):
            p = each + 1
            print(p, exon, sp)
            try:
                alist = stacks[p]
            except:
                newSeq[0].append('-')
                newSeq[1].append('-')
                continue

            for i in range(1):
                chars = [each[i] for each in alist]
                mainChar1, mainChar2 = tell_dominate_char_2phase(chars, cov)
                newSeq[0].append(mainChar1)
                newSeq[1].append(mainChar2)

        for i, eachSeq in enumerate(newSeq):
            line = '%s_%s\t%s\t%s\n' % (exon, i+1, sp, ''.join(eachSeq))
            new.append(line)

cmn.write_file(''.join(new), dn)


