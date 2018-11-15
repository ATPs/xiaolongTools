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

map_dict = {"TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"STOP", "TAG":"STOP",
    "TGT":"C", "TGC":"C", "TGA":"STOP", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",}

def check_translation(dna):
    #take the longest
    #the first element to tell if it is possitive. True = positive; False = reverse
    longest_frame = (True, (0,0))
    for i in range(3):
        frame = check_transcript_frame(dna[i:])
        if frame[1] - frame[0] > (longest_frame[1][1] - longest_frame[1][0]):
            longest_frame = (True, [x+i for x in frame])

    for i in range(3):
        reverse = dna[::-1][i:]
        frame = check_transcript_frame(reverse)
        if frame[1] - frame[0] > (longest_frame[1][1] - longest_frame[1][0]):
            longest_frame = (False, [x+i for x in frame])

    return longest_frame

def check_transcript_frame(seq):
    longest = (0, (0, 0))
    i = 0
    j = 0
    for j in range(0, len(seq), 3):
        codon = seq[j:j+3].upper()
        try:
            char = map_dict[codon]
        except KeyError:
            char = 'X'
        if char == 'STOP':
            length = j - i
            if length > longest[0]:
                longest = (length, (i, j))
            i = j + 3
    if char != 'STOP':
        #no stop codon at last
        length = j - i
        if length > longest[0]:
            longest = (length, (i, j))
    return (i, j)


def dna2transcripts(seq, isSTOP=True):
    new = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i+3].upper()
        try:
            char = map_dict[codon]
        except KeyError:
            char = 'X'
        new.append(char)

    transcripts = []
    #find the longest strength
    i = 0
    j = 0
    current = []
    while(j<len(new)):
        char = new[j]
        #if char != 'STOP' and char != 'X':
        if char != 'STOP':
            current.append(char)
        elif char == 'STOP':
            transcripts.append(''.join(current))
            current = []
        j += 1
    if len(current) != 0:
        transcripts.append(''.join(current))

    return transcripts

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def translate_dna_frame(seq, frame):
    if not frame[0]:
        seq = seq[::-1]

    i, j = frame[1]
    protein = []
    for index in range(i, j+3, 3):
        codon = seq[index:index+3].upper()
        try:
            char = map_dict[codon]
        except KeyError:
            char = 'X'

        #very conservative here
        #don't stop because I don't know if this is a sequencing error
        if char == 'STOP':
            char = 'X'
        protein.append(char)
    return ''.join(protein)




if __name__=='__main__':
    #options=parse_options()
    try:
        fn, ref_sp = sys.argv[1:3]
    except:
        print("Usage: *.py fasta refID", file=sys.stderr)
        sys.exit()

    #firstly translate 3935 to get the reference frame
    #then use this frame to translate others

    #ref_sp = '3935'
    seqDict = read_fa(fn)

    refseqs = [seqDict[name] for name in seqDict
            if name.split('_')[0] == ref_sp]

    frame = check_translation(refseqs[0])

    new = []
    for name in seqDict:
        dna = seqDict[name]
        protein = translate_dna_frame(dna, frame)
        new.append('>%s\n%s\n' % (name, protein))

    print(''.join(new))





