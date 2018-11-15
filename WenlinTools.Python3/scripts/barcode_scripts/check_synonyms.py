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
def parse_paired_fa(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        if line[0] == '>':
            name = line[1:]
            if '(assembled)' in name:
                key = name
                isAss = True
                if key not in adict:
                    adict[key] = [None, {}]
            else:
                isAss = False
        else:
            if isAss:
                adict[key][0] = line
            else:
                adict[key][1][name] = line
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mapDict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
            "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
                "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
                    "TGT":"C", "TGC":"C", "TGA":"W", "TGG":"W",
                        "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
                            "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
                                "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
                                    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
                                        "ATT":"I", "ATC":"I", "ATA":"M", "ATG":"M",
                                            "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
                                                "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
                                                    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
                                                        "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
                                                            "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
                                                                "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
                                                                    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}




def dna2prot(seq):
    new = []
    for i in range(0, len(seq), 3):
        char = seq[i:i+3]
        try:
            prot = mapDict[char]
        except KeyError:
            prot = 'X'

        new.append(prot)
    return ''.join(new)

def compare_proteins(seq1, seq2):
    alist = []
    i = 0
    for char1, char2 in zip(seq1, seq2):
        if char1 == 'X' or char2 == 'X':
            i += 1
            continue

        if char1 != char2:
            alist.append('%s--%s:%s' % (i, char1, char2))
        i += 1
    return alist

def check_dna_diff(pos, seq1, seq2):
    pos = int(pos.split('--')[0])
    i = pos * 3
    j = i + 3
    dna1 = seq1[i:j]
    dna2 = seq2[i:j]
    print(dna1, dna2)
    alist = []
    k = 0
    for char1, char2 in zip(dna1, dna2):
        if char1 != char2:
            alist.append('%s--%s:%s' % (k+i, char1, char2))
        k += 1
    return alist


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *paired.fa", file=sys.stderr)
        sys.exit()


    pair_dict = parse_paired_fa(fn)

    for name in pair_dict:
        assembled_seq, compared_seqs = pair_dict[name]
        base_prot = dna2prot(assembled_seq[1:])
        for seqName in compared_seqs:
            seq = compared_seqs[seqName][1:]
            prot = dna2prot(seq)
            #print prot
            #print base_prot
            diffP = compare_proteins(prot, base_prot)
            if len(diffP) != 0:
                for prot_pos in diffP:
                    dna_pos = check_dna_diff(prot_pos, assembled_seq[1:], seq)
                    print('diff positions: %s vs %s; %s %s' % (name, seqName, prot_pos, dna_pos))
        print('#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#')



