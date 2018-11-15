import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
transfer_dict = {
        'CT': 'R',
        'AG': 'Y',
        'AT': 'W',
        'CG': 'S',
        'GT': 'M',
        'AC': 'K'
}

def collapse_seqs(seq1, seq2):
    seq = []
    for i in range(len(seq1)):
        char1, char2 = [seq1[i], seq2[i]]
        if char1 == char2:
            seq.append(char1)
        else:
            chars = [char1, char2]
            chars.sort()
            char = transfer_dict[''.join(chars)]
            seq.append(char)

    return ''.join(seq)
        

fn = sys.argv[1]
cutoff = int(sys.argv[2])

seqDict = {}
seq = []
with open(fn) as fp:
    for line in fp:
        if line[0] == '>':
            if seq != []:
                seq = ''.join(seq)
                try:
                    seqDict[sp].append(seq)
                except KeyError:
                    seqDict[sp] = [seq]

            sp = line[1:].split('_')[0]
            seq = []
        else:
            seq.append(line.strip())

#last seq
seqDict[sp].append(''.join(seq))

new = []

for sp in seqDict:
    seq1, seq2 = seqDict[sp]
    diffN = sum([seq1[i] != seq2[i]
            for i in range(len(seq1))])

    if diffN < cutoff:
        if diffN == 0:
            seq = seq1
            defline = '%s_unique' % sp
        else:    
            seq = collapse_seqs(seq1, seq2)
            defline = '%s_diff%s' % (sp, diffN)
        fasta = '>%s\n%s\n' % (defline, seq)
    else:            
        #need to keep both copy
        label = '%s_diff%s' % (sp, diffN)
        fasta = '>%s_cp1\n%s\n>%s_cp2\n%s\n' % (label, seq1, label, seq2)
    new.append(fasta)

dn = '%s_collapse_cut%s.fa' % (cmn.lastName(fn).replace('.fa', '') , cutoff)
cmn.write_file(''.join(new), dn)
