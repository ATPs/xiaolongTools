#!/usr/bin/python
import os, sys, numpy, string

python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn


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




try:
    fblast, fread = sys.argv[1:3]
except:
    print('*.py blastExon.dict.pkl readSp.dict.pkl', file=sys.stderr)
    sys.exit()


blast_dict = cmn.pickle_read(fblast)

#reads[sp][name] = seq
reads = cmn.pickle_read(fread)
splist = set(reads.keys())

for exon in blast_dict:
    info = blast_dict[exon]

    #filter by the splist to avoid error
    info = [line for line in info
            if line.split()[0] in splist]

    if len(info) == 0:
        print('missing regions for %s in %s' % (exon, splist), file=sys.stderr)
        continue

    map = {}
    ref = {}
    hitcount = {}
    LENGTH = 0
    for line in info:
        words = line.split()
        sp = words[0]
        rd = words[2]
        qstart = int(words[6])
        if not LENGTH:
            LENGTH = int(words[5])
        qseq = words[11]
        sstart = int(words[9])
        if int(words[10]) > sstart:
            sdir = "+"
        else:
            sdir = "-"
        sseq = words[12].replace('*', '-')
        if len(qseq) - qseq.count("-") >= 20:
            for i in range(len(qseq) - 19):
                qmseg = qseq[i:i+20]
                smseg = sseq[i:i+20]
                if len(qmseg) == 20 and not "-" in qmseg and not "-" in smseg:
                    qbseg = qseq[:i]
                    qbegin = qstart + len(qbseg) - qbseg.count("-")
                    sbseg = sseq[:i]
                    if sdir == "+":
                        sbegin = sstart + 3*(len(sbseg) - sbseg.count("-"))
                        send = sbegin + 59
                    elif sdir == "-":
                        sbegin = sstart - 3*(len(sbseg) - sbseg.count("-"))
                        send = sbegin - 59

                    try:
                        ref[qbegin]
                    except KeyError:
                        ref[qbegin] = qmseg
                    try:
                        map[qbegin]
                    except KeyError:
                        map[qbegin] = {}
                    try:
                        map[qbegin][sp]
                    except KeyError:
                        map[qbegin][sp] = {}
                    try:
                        map[qbegin][sp][rd].append([sbegin,send,smseg])
                    except KeyError:
                        map[qbegin][sp][rd] = [[sbegin,send,smseg]]

                    try:
                        hitcount[qbegin]
                    except KeyError:
                        hitcount[qbegin] = {}
                    try:
                        hitcount[qbegin][smseg] += 1
                    except KeyError:
                        hitcount[qbegin][smseg] = 1

    comp = {}
    comp['A'] = 'T'
    comp['T'] = 'A'
    comp['G'] = 'C'
    comp['C'] = 'G'
    comp['N'] = 'N'

    goodmap = {}
    for qstart in list(hitcount.keys()):
        cutoff = numpy.median(list(hitcount[qstart].values()))
        bestseq = ""
        bestvalue = 20
        for seq in list(hitcount[qstart].keys()):
            if hitcount[qstart][seq] >= cutoff:
                refseq = ref[qstart]
                countdiff = diff_letters(seq,refseq)
                if countdiff < bestvalue:
                    bestvalue = countdiff
                    bestseq = seq
        if bestseq:
            for sp in list(map[qstart].keys()):
                dnaseqs = []
                for rd in list(map[qstart][sp].keys()):
                    for item in map[qstart][sp][rd]:
                        start = item[0]
                        end = item[1]
                        seq = item[2]
                        if diff_letters(seq,bestseq) <= 2:
                            try:
                                goodmap[sp]
                            except KeyError:
                                goodmap[sp] = {}
                            try:
                                goodmap[sp][rd]
                            except KeyError:
                                goodmap[sp][rd] = {}
                            if start < end:
                                dnaseq = reads[sp][rd][start-1:end]
                            else:
                                tmpseq = reads[sp][rd][end-1:start]
                                dnaseq = ""
                                for char in tmpseq[::-1]:
                                    dnaseq += comp[char]

                            for qposi in range(qstart,qstart+20):
                                codon = dnaseq[(qposi-qstart)*3:(qposi-qstart)*3+3]
                                try:
                                    goodmap[sp][rd][qposi]
                                    if codon != goodmap[sp][rd][qposi]:
                                        goodmap[sp][rd][qposi] = "---"
                                except KeyError:
                                    goodmap[sp][rd][qposi] = codon

    seqs = {}
    for i in range(1,LENGTH+1):
        for sp in list(goodmap.keys()):
            try:
                seqs[sp]
            except KeyError:
                seqs[sp] = {}
            for rd in list(goodmap[sp].keys()):
                try:
                    seqs[sp][rd]
                except KeyError:
                    seqs[sp][rd] = ""
                try:
                    seqs[sp][rd] += goodmap[sp][rd][i]
                except KeyError:
                    seqs[sp][rd] += "---"
    
    #print seqs

    cseqs = {}
    for sp in list(seqs.keys()):
        cseqs[sp + "_1"] = ""
        #cseqs[sp + "_2"] = ""
        for i in range(0,LENGTH*3):
            nts = []
            for rd in list(seqs[sp].keys()):
                #print i, sp, rd, seqs[sp][rd]
                #try:
                #    char = seqs[sp][rd][i]
                #except:
                #    continue
                char = seqs[sp][rd][i]
                if char != "-":
                    nts.append(char)
            if not nts:
                cseqs[sp + "_1"] += "-"
                #cseqs[sp + "_2"] += "-"
            else:
                good_nts = most_common(nts)
                #cseqs[sp + "_1"] += good_nts[0]
                cseqs[sp + "_1"] += good_nts
                #cseqs[sp + "_2"] += good_nts[1]
        #print sp + "_1" + "\t" + cseqs[sp + "_1"]
        #print sp + "_2" + "\t" + cseqs[sp + "_2"]

        print('%s\t%s\t%s' % (exon, sp, cseqs[sp + "_1"]))
        #print '%s\t%s_2\t%s' % (exon, sp, cseqs[sp + "_2"])
