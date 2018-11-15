#a lib to keep most of the function for sequencing data
import math
dna_gaps = set(['-', 'N'])

def compute_sequence_entropy(seq, kmer_size=20):
    seq = seq.upper()
    #currently, only deal with dna and don't consider gaps
    for char in dna_gaps:
        seq = seq.replace(char, '')

    length = len(seq)
    count_dict = {}
    for i in range(length - kmer_size):
        seg = seq[i:i+kmer_size]
        try:
            count_dict[seg] += 1
        except KeyError:
            count_dict[seg] = 1


    score = 0
    for key in count_dict:
        count = count_dict[key]
        P = float(count) / (length - kmer_size)
        score += P * math.log(P)
        #score += P * math.log(P, 10)

    return score
