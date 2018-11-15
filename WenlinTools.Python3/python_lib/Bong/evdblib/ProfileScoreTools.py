import sys
import math
import array
import os

##############################################################
# Parse MSA file
##############################################################
#simple MSA reader,
# This function assumes that the MSA contains first line id as
# 'QUERY' and end flanked bla..bla...
# Basically the format for Ruslans 97% sequence purger...


def parse_msa(fp):
    seqs = []
    lc = 0
    for l in fp.readlines():
        if l == '\n':
            continue
        l = l.split()
        #in case of 'QUERY' line
        #reset line counter (lc) to 0
        if l[0] == 'QUERY':
            lc = 0

        if len(seqs) < lc + 1:
            seqs.append(l[-1])
        else:
            seqs[lc] += l[-1]
        lc += 1

    return seqs


#simplest form of MSA input
#multiple lines where each line only contains sequence
#each sequences are in one line
def parse_hua_msa(fp):
    seqs = []
    for l in fp.readlines():
        seqs.append(l[:-1])
    return seqs

#simple score chooser from score list.
#this function is used for selecting a representing score from list of shift random scores
##########################
# Shift selection funtion
##########################


def select_random_score(scores, type):
    if type == 'min':
        return min(scores)
    elif type == 'avg':
        return 1.0 * sum(scores) / len(scores)
    elif type == 'trimmed_avg':
        scores.sort()
        n = max(1, int(len(scores) * 0.9))
        return 1.0 * sum(scores[:n]) / n

    elif type == 'median':
        scores.sort()
        n = int(len(scores) * 0.5)
        if len(scores) % 2:
            return 0.5 * (scores[n] + scores[n - 1])
        else:
            return scores[n]
    else:
        print('Error! shift selection type is not valid', type, file=sys.stderr)
        sys.exit()


#########################
#calculate Neff
#########################
def calculate_neff(msa, pos, aa):
    #selecting submsa's for having same amino acid
    #as given aa in the given column (pos)
    submsa_i = []
    for i, seq in enumerate(msa):
        if seq[pos] == aa:
            submsa_i.append(i)

    if not submsa_i:
        return 0.0
    #if len(submsa_i) == 1 :
    #        return 1.0

    #count of column
    #where any position is non-gap
    col_count = 0
    sum_diff = 0
    for i in range(len(msa[0])):
        #first approach to get number of
        #different types of amino acids
        aa_diff = set()
        for j in submsa_i:
            #if msa[j][i].isalpha() :
            #count only non-gap and 20 standard amino acids
            #if msa[j][i] in aa_list3 :
            aa_diff.add(msa[j][i])
            '''
            if len(aa_diff) == 20 :
                    aa_diff = aa_diff.intersection(aa_list3)
                    if len(aa_diff) == 20 :
                            break
            '''

        aa_diff = aa_diff.intersection(aa_list3)

        if aa_diff:
            sum_diff += len(aa_diff)
            col_count += 1
            #print col_count, aa_diff

    if col_count == 0:
        return 0.0

    #debug
    #print sum_diff

    return math.log(1 - ((0.05 * sum_diff) / col_count)) / math.log(0.95)

##########
# constant matrices
###########
aa_list = "WFYMLIVACGPTSNQDEHRK"
aa_list2 = "-WFYMLIVACGPTSNQDEHRK"
aa_list3 = set("WFYMLIVACGPTSNQDEHRK")
#Robinson background probabilities
p_rbnsn = {'W': 0.013298, 'F': 0.038556, 'Y': 0.032165, 'M': 0.022425, 'L': 0.090191, 'I': 0.05142, 'V': 0.064409, 'A': 0.078047, 'C': 0.019246, 'G': 0.073772, 'P': 0.052028, 'T': 0.058413, 'S': 0.071198, 'N': 0.044873, 'Q': 0.042644, 'D': 0.05364, 'E': 0.062949, 'H': 0.021992, 'R': 0.051295, 'K': 0.057438}

q_blosum62 = {
    'W': {"W": 0.0065, "F": 0.0008, "Y": 0.0009, "M": 0.0002, "L": 0.0007, "I": 0.0004, "V": 0.0004, "A": 0.0004, "C": 0.0001, "G": 0.0004, "P": 0.0001, "T": 0.0003, "S": 0.0003, "N": 0.0002, "Q": 0.0002, "D": 0.0002, "E": 0.0003, "H": 0.0002, "R": 0.0003, "K": 0.0003, },
'F': {"W": 0.0008, "F": 0.0183, "Y": 0.0042, "M": 0.0012, "L": 0.0054, "I": 0.0030, "V": 0.0026, "A": 0.0016, "C": 0.0005, "G": 0.0012, "P": 0.0005, "T": 0.0012, "S": 0.0012, "N": 0.0008, "Q": 0.0005, "D": 0.0008, "E": 0.0009, "H": 0.0008, "R": 0.0009, "K": 0.0009, },
'Y': {"W": 0.0009, "F": 0.0042, "Y": 0.0102, "M": 0.0006, "L": 0.0022, "I": 0.0014, "V": 0.0015, "A": 0.0013, "C": 0.0003, "G": 0.0008, "P": 0.0005, "T": 0.0009, "S": 0.0010, "N": 0.0007, "Q": 0.0007, "D": 0.0006, "E": 0.0009, "H": 0.0015, "R": 0.0009, "K": 0.0010, },
'M': {"W": 0.0002, "F": 0.0012, "Y": 0.0006, "M": 0.0040, "L": 0.0049, "I": 0.0025, "V": 0.0023, "A": 0.0013, "C": 0.0004, "G": 0.0007, "P": 0.0004, "T": 0.0010, "S": 0.0009, "N": 0.0005, "Q": 0.0007, "D": 0.0005, "E": 0.0007, "H": 0.0004, "R": 0.0008, "K": 0.0009, },
'L': {"W": 0.0007, "F": 0.0054, "Y": 0.0022, "M": 0.0049, "L": 0.0371, "I": 0.0114, "V": 0.0095, "A": 0.0044, "C": 0.0016, "G": 0.0021, "P": 0.0014, "T": 0.0033, "S": 0.0024, "N": 0.0014, "Q": 0.0016, "D": 0.0015, "E": 0.0020, "H": 0.0010, "R": 0.0024, "K": 0.0025, },
'I': {"W": 0.0004, "F": 0.0030, "Y": 0.0014, "M": 0.0025, "L": 0.0114, "I": 0.0184, "V": 0.0120, "A": 0.0032, "C": 0.0011, "G": 0.0014, "P": 0.0010, "T": 0.0027, "S": 0.0017, "N": 0.0010, "Q": 0.0009, "D": 0.0012, "E": 0.0012, "H": 0.0006, "R": 0.0012, "K": 0.0016, },
'V': {"W": 0.0004, "F": 0.0026, "Y": 0.0015, "M": 0.0023, "L": 0.0095, "I": 0.0120, "V": 0.0196, "A": 0.0051, "C": 0.0014, "G": 0.0018, "P": 0.0012, "T": 0.0036, "S": 0.0024, "N": 0.0012, "Q": 0.0012, "D": 0.0013, "E": 0.0017, "H": 0.0006, "R": 0.0016, "K": 0.0019, },
'A': {"W": 0.0004, "F": 0.0016, "Y": 0.0013, "M": 0.0013, "L": 0.0044, "I": 0.0032, "V": 0.0051, "A": 0.0215, "C": 0.0016, "G": 0.0058, "P": 0.0022, "T": 0.0037, "S": 0.0063, "N": 0.0019, "Q": 0.0019, "D": 0.0022, "E": 0.0030, "H": 0.0011, "R": 0.0023, "K": 0.0033, },
'C': {"W": 0.0001, "F": 0.0005, "Y": 0.0003, "M": 0.0004, "L": 0.0016, "I": 0.0011, "V": 0.0014, "A": 0.0016, "C": 0.0119, "G": 0.0008, "P": 0.0004, "T": 0.0009, "S": 0.0010, "N": 0.0004, "Q": 0.0003, "D": 0.0004, "E": 0.0004, "H": 0.0002, "R": 0.0004, "K": 0.0005, },
'G': {"W": 0.0004, "F": 0.0012, "Y": 0.0008, "M": 0.0007, "L": 0.0021, "I": 0.0014, "V": 0.0018, "A": 0.0058, "C": 0.0008, "G": 0.0378, "P": 0.0014, "T": 0.0022, "S": 0.0038, "N": 0.0029, "Q": 0.0014, "D": 0.0025, "E": 0.0019, "H": 0.0010, "R": 0.0017, "K": 0.0025, },
'P': {"W": 0.0001, "F": 0.0005, "Y": 0.0005, "M": 0.0004, "L": 0.0014, "I": 0.0010, "V": 0.0012, "A": 0.0022, "C": 0.0004, "G": 0.0014, "P": 0.0191, "T": 0.0014, "S": 0.0017, "N": 0.0009, "Q": 0.0008, "D": 0.0012, "E": 0.0014, "H": 0.0005, "R": 0.0010, "K": 0.0016, },
'T': {"W": 0.0003, "F": 0.0012, "Y": 0.0009, "M": 0.0010, "L": 0.0033, "I": 0.0027, "V": 0.0036, "A": 0.0037, "C": 0.0009, "G": 0.0022, "P": 0.0014, "T": 0.0125, "S": 0.0047, "N": 0.0022, "Q": 0.0014, "D": 0.0019, "E": 0.0020, "H": 0.0007, "R": 0.0018, "K": 0.0023, },
'S': {"W": 0.0003, "F": 0.0012, "Y": 0.0010, "M": 0.0009, "L": 0.0024, "I": 0.0017, "V": 0.0024, "A": 0.0063, "C": 0.0010, "G": 0.0038, "P": 0.0017, "T": 0.0047, "S": 0.0126, "N": 0.0031, "Q": 0.0019, "D": 0.0028, "E": 0.0030, "H": 0.0011, "R": 0.0023, "K": 0.0031, },
'N': {"W": 0.0002, "F": 0.0008, "Y": 0.0007, "M": 0.0005, "L": 0.0014, "I": 0.0010, "V": 0.0012, "A": 0.0019, "C": 0.0004, "G": 0.0029, "P": 0.0009, "T": 0.0022, "S": 0.0031, "N": 0.0141, "Q": 0.0015, "D": 0.0037, "E": 0.0022, "H": 0.0014, "R": 0.0020, "K": 0.0024, },
'Q': {"W": 0.0002, "F": 0.0005, "Y": 0.0007, "M": 0.0007, "L": 0.0016, "I": 0.0009, "V": 0.0012, "A": 0.0019, "C": 0.0003, "G": 0.0014, "P": 0.0008, "T": 0.0014, "S": 0.0019, "N": 0.0015, "Q": 0.0073, "D": 0.0016, "E": 0.0035, "H": 0.0010, "R": 0.0025, "K": 0.0031, },
'D': {"W": 0.0002, "F": 0.0008, "Y": 0.0006, "M": 0.0005, "L": 0.0015, "I": 0.0012, "V": 0.0013, "A": 0.0022, "C": 0.0004, "G": 0.0025, "P": 0.0012, "T": 0.0019, "S": 0.0028, "N": 0.0037, "Q": 0.0016, "D": 0.0213, "E": 0.0049, "H": 0.0010, "R": 0.0016, "K": 0.0024, },
'E': {"W": 0.0003, "F": 0.0009, "Y": 0.0009, "M": 0.0007, "L": 0.0020, "I": 0.0012, "V": 0.0017, "A": 0.0030, "C": 0.0004, "G": 0.0019, "P": 0.0014, "T": 0.0020, "S": 0.0030, "N": 0.0022, "Q": 0.0035, "D": 0.0049, "E": 0.0161, "H": 0.0014, "R": 0.0027, "K": 0.0041, },
'H': {"W": 0.0002, "F": 0.0008, "Y": 0.0015, "M": 0.0004, "L": 0.0010, "I": 0.0006, "V": 0.0006, "A": 0.0011, "C": 0.0002, "G": 0.0010, "P": 0.0005, "T": 0.0007, "S": 0.0011, "N": 0.0014, "Q": 0.0010, "D": 0.0010, "E": 0.0014, "H": 0.0093, "R": 0.0012, "K": 0.0012, },
'R': {"W": 0.0003, "F": 0.0009, "Y": 0.0009, "M": 0.0008, "L": 0.0024, "I": 0.0012, "V": 0.0016, "A": 0.0023, "C": 0.0004, "G": 0.0017, "P": 0.0010, "T": 0.0018, "S": 0.0023, "N": 0.0020, "Q": 0.0025, "D": 0.0016, "E": 0.0027, "H": 0.0012, "R": 0.0178, "K": 0.0062, },
'K': {"W": 0.0003, "F": 0.0009, "Y": 0.0010, "M": 0.0009, "L": 0.0025, "I": 0.0016, "V": 0.0019, "A": 0.0033, "C": 0.0005, "G": 0.0025, "P": 0.0016, "T": 0.0023, "S": 0.0031, "N": 0.0024, "Q": 0.0031, "D": 0.0024, "E": 0.0041, "H": 0.0012, "R": 0.0062, "K": 0.0161, }
}


def calculate_mean_number_of_symbols(msa):
    #mean number of different symbols in the given position
    Nc = 0.0
    line_count = 0
    for i in range(len(msa[0])):
        aa_set = set()
        for j in range(len(msa)):
            aa_set.add(msa[j][i])
        #aa_set.discard( '-' )
        #'''
        if aa_set:
            Nc += len(aa_set)
            line_count += 1
        #'''
        '''
        if len(aa_set)>1 :
                line_count += 1
                Nc += len(aa_set)
        elif len(aa_set) == 1 :
                if not '-' in aa_set :
                        line_count += 1
                        Nc += 1
        '''

    Nc /= line_count
    return Nc


def calculate_mean_number_of_symbols2(neffs):
    col_count = 0
    symbol_count = 0
    for effs in neffs:
        col_flag = 0
        for a in aa_list2:
            if effs[a] > 0.0:
                symbol_count += 1
                col_flag = 1
        if col_flag:
            col_count += 1

    return (1.0 * symbol_count) / col_count


def calculate_effective_number_of_sequences(neffs):
    '''
    This function calculates effective number of sequences
    for a given Neffs of a MSA.

    This is implemented by a simple avg of summation of
    Neffs for each column.
    '''
    l = [sum(neff.values()) for neff in neffs]
    return sum(l) / len(l)

##########
# target frequency
##########
# f is effective frequency
# pos is the column position
# msa is the MSA sequences
# aa is amino acid type


def calculate_target_frequency(msa, feffs, pos, aa, alpha):
    #calculate pseudocount
    g = 0.0
    for a in aa_list:
        g += feffs[pos][a] * q_blosum62[a][aa] / p_rbnsn[a]

    return (alpha * feffs[pos][aa] + 10 * g) / (alpha + 10)

############
# pssm
############
#def calculate_pssm( Q, p ) :
    #print 'pssm', Q, p,
    #print math.log( Q/p )
#       return math.log( Q/p )

####################################################
# profile scoring schemes
####################################################

##########
# COMPASS
##########


def cal_compass_like_score(eq1, eq2, neffs1, neffs2, pssm1, pssm2):
    score = 0.0
    for i, j in zip(eq1, eq2):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting compass score 0 for the column
        if not pssm1[i] or not pssm2[j]:
            continue

        sum_neff1 = 0.0
        sum_neff2 = 0.0
        for a in aa_list:
            sum_neff1 += neffs1[i][a]
            sum_neff2 += neffs2[j][a]

        if sum_neff1 + sum_neff2 == 2:
            c1 = 1.0
            c2 = 1.0
        else:
            c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
            c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

        s1 = 0.0
        s2 = 0.0
        for a in aa_list:
            s1 += neffs1[i][a] * pssm2[j][a]
            s2 += neffs2[j][a] * pssm1[i][a]

        score += c1 * s1 + c2 * s2

    return score


def cal_compass_like_shuffle_random_score(eq1, eq2, neffs1, neffs2, pssm1, pssm2):
    score = 0.0
    count = 0
    for i in eq1:
        for j in eq2:
            #if one of the column do not have sequences
            # it should be skipped.. -> effectively getting compass score 0 for the column
            if not pssm1[i] or not pssm2[j]:
                continue
            count += 1

            sum_neff1 = 0.0
            sum_neff2 = 0.0
            for a in aa_list:
                sum_neff1 += neffs1[i][a]
                sum_neff2 += neffs2[j][a]

            if sum_neff1 + sum_neff2 == 2:
                c1 = 1.0
                c2 = 1.0
            else:
                c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
                c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

            s1 = 0.0
            s2 = 0.0
            for a in aa_list:
                s1 += neffs1[i][a] * pssm2[j][a]
                s2 += neffs2[j][a] * pssm1[i][a]

            score += c1 * s1 + c2 * s2

    return score / math.sqrt(count)

##########
# COMPRASS LIKE SCORE
# Not yet implemented !!!
##########


def cal_comprass_like_score(eq1, eq2, sec1, sec2, neffs1, neffs2, pssm1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting compass score 0 for the column
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue

        sum_neff1 = 0.0
        sum_neff2 = 0.0
        for a in aa_list:
            sum_neff1 += neffs1[eq1[i]]
            sum_neff2 += neffs2[eq2[i]]
        c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
        c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

        s1 = 0.0
        s2 = 0.0
        for a in aa_list:
            s1 += neffs1[eq1[i]][a] * pssm2[eq2[i]][a]
            s2 += neffs2[eq2[i]][a] * pssm1[eq1[i]][a]

        score += c1 * s1 + c2 * s2

    return score


##########
# Pearson's Correlation
##########
def cal_pearson_correlation(eq1, eq2, pssm1, pssm2):
    score = 0.0

    for i in range(len(eq1)):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting pearson score 0 for the column
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue

        #get mean number for each columns
        wbar1 = 0.0
        wbar2 = 0.0
        for a in aa_list:
            wbar1 += pssm1[eq1[i]][a]
            wbar2 += pssm2[eq2[i]][a]
        wbar1 = wbar1 / 20.0
        wbar2 = wbar2 / 20.0

        p12 = 0.0
        p11 = 0.0
        p22 = 0.0

        for a in aa_list:
            p12 += (pssm1[eq1[i]][a] - wbar1) * (pssm2[eq2[i]][a] - wbar2)
            p11 += (pssm1[eq1[i]][a] - wbar1) * (pssm1[eq1[i]][a] - wbar1)
            p22 += (pssm2[eq2[i]][a] - wbar2) * (pssm2[eq2[i]][a] - wbar2)

        score += p12 / math.sqrt(p11 * p22)

    return score


##########
# PICASSO3Q
##########
def cal_picasso3q(eq1, eq2, Qs1, Qs2, pssm1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting picasso score 0 for the column
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue

        s1 = 0.0
        s2 = 0.0
        for a in aa_list:
            s1 += Qs1[eq1[i]][a] * pssm2[eq2[i]][a]
            s2 += Qs2[eq2[i]][a] * pssm1[eq1[i]][a]

        score += s1 + s2

    return score


##########
# Cross Product (Sum of pairs, LogAverage)
##########
def cal_cross_product(eq1, eq2, Qs1, Qs2):
    score = 0.0
    for i in range(len(eq1)):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting crossproduct score 0 for the column
        if not Qs1[eq1[i]] or not Qs2[eq2[i]]:
            continue

        s12 = 0.0
        for a in aa_list:
            for b in aa_list:
                s12 += Qs1[eq1[i]][a] * Qs2[eq2[i]][b] * q_blosum62[a][b]
        score += math.log(s12)
    return score


##########
# Dot product (DotPOdds)
##########
def cal_dot_product(eq1, eq2, pssm1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        #if one of the column do not have sequences
        # it should be skipped.. -> effectively getting dot product score 0 for the column
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue

        s12 = 0.0
        for a in aa_list:
            s12 += pssm1[eq1[i]][a] * pssm2[eq2[i]][a]

        score += s12

    return score

##########
# RPS-BLAST Like
##########


def cal_rps_blast_like_score(eq1, eq2, seq1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        if not pssm2[eq2[i]]:
            continue
        score += pssm2[eq2[i]][seq1[eq1[i]]]

    return score

##########
# PSI-BLAST Like
##########


def cal_psi_blast_like_score(eq1, eq2, pssm1, seq2):
    score = 0.0
    for i in range(len(eq1)):
        if not pssm1[eq1[i]]:
            continue
        score += pssm1[eq1[i]][seq2[eq2[i]]]

    return score

##########
# simple combined_blast like
##########


def cal_simple_combined_blast_like_score(eq1, eq2, seq1, seq2, pssm1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue

        #in case of good amino acid character
        if seq1[eq1[i]] in aa_list:
            score += pssm2[eq2[i]][seq1[eq1[i]]]
        if seq2[eq2[i]] in aa_list:
            score += pssm1[eq1[i]][seq2[eq2[i]]]
    return score

##########
# combined_blast like
##########


def cal_combined_blast_like_score(eq1, eq2, seq1, seq2, neffs1, neffs2, pssm1, pssm2):
    score = 0.0
    for i in range(len(eq1)):
        if not pssm1[eq1[i]] or not pssm2[eq2[i]]:
            continue
        sum_neff1 = 0.0
        sum_neff2 = 0.0
        for a in aa_list:
            sum_neff1 += neffs1[eq1[i]][a]
            sum_neff2 += neffs2[eq2[i]][a]
        if sum_neff1 + sum_neff2 == 2:
            c1 = c2 = 1.0
        else:
            c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
            c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

        if seq1[eq1[i]] in aa_list:
            score += c1 * pssm2[eq2[i]][seq1[eq1[i]]]
        if seq2[eq2[i]] in aa_list:
            score += c2 * pssm1[eq1[i]][seq2[eq2[i]]]

    return score

#################
# Faster calculation code by calculate everything together
##################


def cal_all_profile_scores(eq1, eq2, seq1, seq2, neffs1, neffs2, Qs1, Qs2, pssm1, pssm2):

    #premake the sequences and arrays containing only equivalent portions
    temp = []
    for i in range(len(eq1)):
        temp.append(seq1[eq1[i]])
    seq1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(seq2[eq2[i]])
    seq2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(neffs1[eq1[i]])
    neffs1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(neffs2[eq2[i]])
    neffs2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(Qs1[eq1[i]])
    Qs1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(Qs2[eq2[i]])
    Qs2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(pssm1[eq1[i]])
    pssm1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(pssm2[eq2[i]])
    pssm2 = temp

    cblscore = 0.0  # combined blast like score
    blscore = 0.0  # blast like score
    dpscore = 0.0  # dot poduct
    cpscore = 0.0  # cross product
    csscore = 0.0  # compass like score
    psscore = 0.0  # picasso score
    pnscore = 0.0  # pearson correlation coefficient
    for i in range(len(eq1)):
        if not pssm1[i] or not pssm2[i]:
            continue

        #processing portion for combined-blast-like score
        sum_neff1 = 0.0
        sum_neff2 = 0.0
        for a in aa_list:
            sum_neff1 += neffs1[i][a]
            sum_neff2 += neffs2[i][a]
        if sum_neff1 + sum_neff2 == 2:
            c1 = c2 = 1.0
        else:
            c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
            c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

        cps12 = 0.0
        for a in aa_list:
            dpscore += pssm1[i][a] * pssm2[i][a]  # dot product
            for b in aa_list:
                cps12 += Qs1[i][a] * Qs2[i][b] * q_blosum62[a][b]
        cpscore += math.log(cps12)

        #blast like scores
        if seq1[i] in aa_list:
            cblscore += c1 * pssm2[i][seq1[i]]
            blscore += pssm2[i][seq1[i]]
        if seq2[i] in aa_list:
            cblscore += c2 * pssm1[i][seq2[i]]
            blscore += pssm1[i][seq2[i]]

        cs1 = 0.0
        cs2 = 0.0
        for a in aa_list:
            cs1 += neffs1[i][a] * pssm2[i][a]  # compass
            cs2 += neffs2[i][a] * pssm1[i][a]
            psscore += Qs1[i][a] * pssm2[i][a]  # picasso
            psscore += Qs2[i][a] * pssm1[i][a]

        csscore += c1 * cs1 + c2 * cs2

        #get mean number for each columns
        wbar1 = 0.0
        wbar2 = 0.0
        for a in aa_list:
            wbar1 += pssm1[i][a]
            wbar2 += pssm2[i][a]
        wbar1 = wbar1 / 20.0
        wbar2 = wbar2 / 20.0

        p12 = 0.0
        p11 = 0.0
        p22 = 0.0

        for a in aa_list:
            p12 += (pssm1[i][a] - wbar1) * (pssm2[i][a] - wbar2)
            p11 += (pssm1[i][a] - wbar1) * (pssm1[i][a] - wbar1)
            p22 += (pssm2[i][a] - wbar2) * (pssm2[i][a] - wbar2)

        pnscore += p12 / math.sqrt(p11 * p22)

    return csscore, pnscore, psscore, dpscore, cpscore, blscore, cblscore

#################
# Faster calculation code by calculate everything together
##################


def cal_all_shuffle_random_profile_scores(eq1, eq2, seq1, seq2, neffs1, neffs2, Qs1, Qs2, pssm1, pssm2):

    #premake the sequences and arrays containing only equivalent portions
    temp = []
    for i in range(len(eq1)):
        temp.append(seq1[eq1[i]])
    seq1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(seq2[eq2[i]])
    seq2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(neffs1[eq1[i]])
    neffs1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(neffs2[eq2[i]])
    neffs2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(Qs1[eq1[i]])
    Qs1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(Qs2[eq2[i]])
    Qs2 = temp

    temp = []
    for i in range(len(eq1)):
        temp.append(pssm1[eq1[i]])
    pssm1 = temp

    temp = []
    for i in range(len(eq2)):
        temp.append(pssm2[eq2[i]])
    pssm2 = temp

    cblscore = 0.0  # combined blast like score
    blscore = 0.0  # blast like score
    dpscore = 0.0  # dot poduct
    cpscore = 0.0  # cross product
    csscore = 0.0  # compass like score
    psscore = 0.0  # picasso score
    pnscore = 0.0  # pearson correlation coefficient
    count = 0
    for i in range(len(eq1)):
        for j in range(len(eq2)):
            if not pssm1[i] or not pssm2[j]:
                continue

            count += 1
            #processing portion for combined-blast-like score
            sum_neff1 = 0.0
            sum_neff2 = 0.0
            for a in aa_list:
                sum_neff1 += neffs1[i][a]
                sum_neff2 += neffs2[j][a]
            if sum_neff1 + sum_neff2 == 2:
                c1 = c2 = 1.0
            else:
                c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
                c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

            cps12 = 0.0
            for a in aa_list:
                dpscore += pssm1[i][a] * pssm2[j][a]  # dot product
                for b in aa_list:
                    cps12 += Qs1[i][a] * Qs2[j][b] * q_blosum62[a][b]
            cpscore += math.log(cps12)

            #blast like scores
            if seq1[i] in aa_list:
                cblscore += c1 * pssm2[j][seq1[i]]
                blscore += pssm2[j][seq1[i]]
            if seq2[j] in aa_list:
                cblscore += c2 * pssm1[i][seq2[j]]
                blscore += pssm1[i][seq2[j]]

            cs1 = 0.0
            cs2 = 0.0
            for a in aa_list:
                cs1 += neffs1[i][a] * pssm2[j][a]  # compass
                cs2 += neffs2[j][a] * pssm1[i][a]
                psscore += Qs1[i][a] * pssm2[j][a]  # picasso
                psscore += Qs2[j][a] * pssm1[i][a]

            csscore += c1 * cs1 + c2 * cs2

            #get mean number for each columns
            wbar1 = 0.0
            wbar2 = 0.0
            for a in aa_list:
                wbar1 += pssm1[i][a]
                wbar2 += pssm2[j][a]
            wbar1 = wbar1 / 20.0
            wbar2 = wbar2 / 20.0

            p12 = 0.0
            p11 = 0.0
            p22 = 0.0

            for a in aa_list:
                p12 += (pssm1[i][a] - wbar1) * (pssm2[j][a] - wbar2)
                p11 += (pssm1[i][a] - wbar1) * (pssm1[i][a] - wbar1)
                p22 += (pssm2[j][a] - wbar2) * (pssm2[j][a] - wbar2)

            pnscore += p12 / math.sqrt(p11 * p22)

    if count:
        rc = 1.0 / math.sqrt(count)
    else:
        rc = 0.0
    return csscore * rc, pnscore * rc, psscore * rc, dpscore * rc, cpscore * rc, blscore * rc, cblscore * rc


def cal_all_shuffle_random2_profile_scores(eq1, eq2, seq1, seq2, neffs1, neffs2, Qs1, Qs2, pssm1, pssm2):

    """
    #premake the sequences and arrays containing only equivalent portions
    temp = []
    for i in xrange( len(eq1) ) :
            temp.append( seq1[eq1[i]] )
    seq1 = temp

    temp = []
    for i in xrange( len(eq2) ) :
            temp.append( seq2[eq2[i]] )
    seq2 = temp

    temp = []
    for i in xrange( len(eq1) ) :
            temp.append( neffs1[eq1[i]] )
    neffs1 = temp

    temp = []
    for i in xrange( len(eq2) ) :
            temp.append( neffs2[eq2[i]] )
    neffs2 = temp

    temp = []
    for i in xrange( len(eq1) ) :
            temp.append( Qs1[eq1[i]] )
    Qs1 = temp

    temp = []
    for i in xrange( len(eq2) ) :
            temp.append( Qs2[eq2[i]] )
    Qs2 = temp

    temp = []
    for i in xrange( len(eq1) ) :
            temp.append( pssm1[eq1[i]] )
    pssm1 = temp

    temp = []
    for i in xrange( len(eq2) ) :
            temp.append( pssm2[eq2[i]] )
    pssm2 = temp
    """

    cblscore = 0.0  # combined blast like score
    blscore = 0.0  # blast like score
    dpscore = 0.0  # dot poduct
    cpscore = 0.0  # cross product
    csscore = 0.0  # compass like score
    psscore = 0.0  # picasso score
    pnscore = 0.0  # pearson correlation coefficient
    count = 0
    for i in range(len(eq1)):
        for j in range(len(eq2)):
            if not pssm1[i] or not pssm2[j]:
                continue

            count += 1
            #processing portion for combined-blast-like score
            sum_neff1 = 0.0
            sum_neff2 = 0.0
            for a in aa_list:
                sum_neff1 += neffs1[i][a]
                sum_neff2 += neffs2[j][a]
            if sum_neff1 + sum_neff2 == 2:
                c1 = c2 = 1.0
            else:
                c1 = (sum_neff2 - 1) / (sum_neff1 + sum_neff2 - 2)
                c2 = (sum_neff1 - 1) / (sum_neff1 + sum_neff2 - 2)

            cps12 = 0.0
            for a in aa_list:
                dpscore += pssm1[i][a] * pssm2[j][a]  # dot product
                for b in aa_list:
                    cps12 += Qs1[i][a] * Qs2[j][b] * q_blosum62[a][b]
            cpscore += math.log(cps12)

            #blast like scores
            if seq1[i] in aa_list:
                cblscore += c1 * pssm2[j][seq1[i]]
                blscore += pssm2[j][seq1[i]]
            if seq2[j] in aa_list:
                cblscore += c2 * pssm1[i][seq2[j]]
                blscore += pssm1[i][seq2[j]]

            cs1 = 0.0
            cs2 = 0.0
            for a in aa_list:
                cs1 += neffs1[i][a] * pssm2[j][a]  # compass
                cs2 += neffs2[j][a] * pssm1[i][a]
                psscore += Qs1[i][a] * pssm2[j][a]  # picasso
                psscore += Qs2[j][a] * pssm1[i][a]

            csscore += c1 * cs1 + c2 * cs2

            #get mean number for each columns
            wbar1 = 0.0
            wbar2 = 0.0
            for a in aa_list:
                wbar1 += pssm1[i][a]
                wbar2 += pssm2[j][a]
            wbar1 = wbar1 / 20.0
            wbar2 = wbar2 / 20.0

            p12 = 0.0
            p11 = 0.0
            p22 = 0.0

            for a in aa_list:
                p12 += (pssm1[i][a] - wbar1) * (pssm2[j][a] - wbar2)
                p11 += (pssm1[i][a] - wbar1) * (pssm1[i][a] - wbar1)
                p22 += (pssm2[j][a] - wbar2) * (pssm2[j][a] - wbar2)

            pnscore += p12 / math.sqrt(p11 * p22)

    #rc = 1.0/math.sqrt(count)
    rc = 1.0 / count * len(eq1)
    return csscore * rc, pnscore * rc, psscore * rc, dpscore * rc, cpscore * rc, blscore * rc, cblscore * rc


######################
######################
#
# Profile based LOOP scores
#
######################
######################

#little utility function
def get_mean(neffs, aascore, seq):
    aastring = 'ACDEFGHIKLMNPQRSTVWY'
    if not neffs:
        return 0.0

    sum_neffs = 0.0
    sum_score = 0.0
    for i, col in enumerate(neffs):
        if not col:
            sum_neffs += 1
            sum_score += aascore[seq[i]]
            continue

        for aa in aastring:
            sum_score += col[aa] * aascore[aa]
            sum_neffs += col[aa]

    #correction for zero sum_neffs score
    if sum_neffs == 0.0:
        return 0.0

    return sum_score / sum_neffs

#D Net charge (Klein et al., 1984)
netcharge = {'A': 0.0, 'C': 0.0, 'E': -1.0, 'D': -1.0, 'G': 0.0, 'F': 0.0, 'I': 0.0, 'H': 0.0, 'K': 1.0, 'M': 0.0, 'L': 0.0, 'N': 0.0, 'Q': 0.0, 'P': 0.0, 'S': 0.0, 'R': 1.0, 'T': 0.0, 'W': 0.0, 'V': 0.0, 'Y': 0.0}
#Kyte-Dolittle scale
hydrophobicity_scale = {'A': 1.8, 'C': 2.5, 'E': -3.5, 'D': -3.5, 'G': -0.40000000000000002, 'F': 2.7999999999999998, 'I': 4.5, 'H': -3.2000000000000002, 'K': -3.8999999999999999, 'M': 1.8999999999999999, 'L': 3.7999999999999998, 'N': -3.5, 'Q': -3.5, 'P': -1.6000000000000001, 'S': -0.80000000000000004, 'R': -4.5, 'T': -0.69999999999999996, 'W': -0.90000000000000002, 'V': 4.2000000000000002, 'Y': -1.3}
#D Normalized van der Waals volume (Fauchere et al., 1988)
van_der_waals = {'A': 1.0, 'C': 2.4300000000000002, 'E': 3.7799999999999998, 'D': 2.7799999999999998, 'G': 0.0, 'F': 5.8899999999999997, 'I': 4.0, 'H': 4.6600000000000001, 'K': 4.7699999999999996, 'M': 4.4299999999999997, 'L': 4.0, 'N': 2.9500000000000002, 'Q': 3.9500000000000002, 'P': 2.7200000000000002, 'S': 1.6000000000000001, 'R': 6.1299999999999999, 'T': 2.6000000000000001, 'W': 8.0800000000000001, 'V': 3.0, 'Y': 6.4699999999999998}
#D Normalized frequency of alpha-helix (Chou-Fasman, 1978b)
helix = {'A': 1.4199999999999999, 'C': 0.69999999999999996, 'E': 1.51, 'D': 1.01, 'G': 0.56999999999999995, 'F': 1.1299999999999999, 'I': 1.0800000000000001, 'H': 1.0, 'K': 1.1599999999999999, 'M': 1.45, 'L': 1.21, 'N': 0.67000000000000004, 'Q': 1.1100000000000001, 'P': 0.56999999999999995, 'S': 0.77000000000000002, 'R': 0.97999999999999998, 'T': 0.82999999999999996, 'W': 1.0800000000000001, 'V': 1.0600000000000001, 'Y': 0.68999999999999995}
#D Normalized frequency of beta-sheet (Chou-Fasman, 1978b)
strand = {'A': 0.82999999999999996, 'C': 1.1899999999999999, 'E': 0.37, 'D': 0.54000000000000004, 'G': 0.75, 'F': 1.3799999999999999, 'I': 1.6000000000000001, 'H': 0.87, 'K': 0.73999999999999999, 'M': 1.05, 'L': 1.3, 'N': 0.89000000000000001, 'Q': 1.1000000000000001, 'P': 0.55000000000000004, 'S': 0.75, 'R': 0.93000000000000005, 'T': 1.1899999999999999, 'W': 1.3700000000000001, 'V': 1.7, 'Y': 1.47}
#D Normalized frequency of coil (Tanaka-Scheraga, 1977)
coil = {'A': 0.94499999999999995, 'C': 0.93200000000000005, 'E': 1.014, 'D': 1.3149999999999999, 'G': 2.355, 'F': 0.622, 'I': 0.67300000000000004, 'H': 0.52500000000000002, 'K': 0.94699999999999995, 'M': 1.028, 'L': 0.75800000000000001, 'N': 1.202, 'Q': 0.70399999999999996, 'P': 0.57899999999999996, 'S': 1.1399999999999999, 'R': 0.36399999999999999, 'T': 0.86299999999999999, 'W': 0.77700000000000002, 'V': 0.56100000000000005, 'Y': 0.90700000000000003}
#D Average flexibility indices (Bhaskaran-Ponnuswamy, 1988)
flexibility_index = {'A': 0.35699999999999998, 'C': 0.34599999999999997, 'E': 0.497, 'D': 0.51100000000000001, 'G': 0.54400000000000004, 'F': 0.314, 'I': 0.46200000000000002, 'H': 0.32300000000000001, 'K': 0.46600000000000003, 'M': 0.29499999999999998, 'L': 0.36499999999999999, 'N': 0.46300000000000002, 'Q': 0.49299999999999999, 'P': 0.50900000000000001, 'S': 0.50700000000000001, 'R': 0.52900000000000003, 'T': 0.44400000000000001, 'W': 0.30499999999999999, 'V': 0.38600000000000001, 'Y': 0.41999999999999998}

#hydrophobicity score


def hydrophobicity(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, hydrophobicity_scale, seq1)
    mean2 = get_mean(neffs2, hydrophobicity_scale, seq2)

    return abs(mean1 - mean2)

#amino acide size loop score


def size(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, van_der_waals, seq1)
    mean2 = get_mean(neffs2, van_der_waals, seq2)

    return abs(mean1 - mean2)

#charge loop score


def charge(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, netcharge, seq1)
    mean2 = get_mean(neffs2, netcharge, seq2)

    return abs(mean1 - mean2)

#length loop socre


def looplength1(seq1, seq2, neffs1, neffs2):
    return abs(len(seq1) - len(seq2))


def looplength2(seq1, seq2, neffs1, neffs2):
    len1 = 0
    for i in neffs1:
        if not i:
            continue
        if 1.0 * i['-'] / sum(i.values()) > 0.5:
            continue
        len1 += 1

    len2 = 0
    for i in neffs2:
        if not i:
            continue
        if 1.0 * i['-'] / sum(i.values()) > 0.5:
            continue
        len2 += 1

    return abs(len1 - len2)

#propensity scores


def helixpropensity(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, helix, seq1)
    mean2 = get_mean(neffs2, helix, seq2)

    return abs(mean1 - mean2)


def coilpropensity(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, coil, seq1)
    mean2 = get_mean(neffs2, coil, seq2)

    return abs(mean1 - mean2)


def strandpropensity(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, strand, seq1)
    mean2 = get_mean(neffs2, strand, seq2)

    return abs(mean1 - mean2)


#flexibility scores

def flexibility(seq1, seq2, neffs1, neffs2):
    mean1 = get_mean(neffs1, flexibility_index, seq1)
    mean2 = get_mean(neffs2, flexibility_index, seq2)

    return abs(mean1 - mean2)


###########################
# loop score interface
###########################
#loop score calculation
#it has correction random score calculation for nloop is 1
def print_profile_based_loop_scores(seq1, seq2, neffs1, neffs2, ploops, cloops, fp):
    score_functions = [hydrophobicity, flexibility, strandpropensity, coilpropensity, helixpropensity, looplength2, looplength1, charge, size]

    for fn in score_functions:
        scores = []
        nloops = len(ploops)
        for i in range(nloops):
            scores.append([])
            for j in range(nloops):
                scores[i].append(fn(seq1[ploops[i][0]:ploops[i][1]], seq2[cloops[j][0]:cloops[j][1]], neffs1[ploops[i][0]:ploops[i][1]], neffs2[cloops[j][0]:cloops[j][1]]))

        scores_r = [0.0] * nloops

        for i in range(nloops):
            for j in range(nloops):
                if i != j:
                    scores_r[i] += scores[i][j]
                    scores_r[j] += scores[i][j]

        # correction for no loop case :
        if not nloops:
            print(fn.__name__ + 'a:', 1.0, 1.0, file=fp)
            print(fn.__name__ + 'b:', 0.0, 0.0, 0.0, 1.0, file=fp)
            continue

        ######################
        # Very important correction for
        # random loop score when there are only one loops
        # in the alignment
        ######################
        if nloops == 1:
            scores_r[0] = (fn(seq1[ploops[0][0]:ploops[0][1]], '', neffs1[ploops[0][0]:ploops[0][1]], []) + fn(seq2[cloops[0][0]:cloops[0][1]], '', neffs2[cloops[0][0]:cloops[0][1]], [])) * 0.5

        else:
            for i in range(nloops):
                scores_r[i] = scores_r[i] / (2.0 * nloops - 2)

        sum_score = 0.0
        sum_score_r = 0.0
        for i in range(nloops):
            sum_score += scores[i][i]
            sum_score_r += scores_r[i]

        scorea = 0.0
        for i in range(nloops):
            try:
                scorea += (scores[i][i] - scores_r[i]) / (0.0 - scores_r[i])
            except ZeroDivisionError:
                pass
        scorea = scorea / nloops

        print(fn.__name__ + 'a:', scorea, file=fp)
        print(fn.__name__ + 'b:', sum_score, 0.0, 0.0, sum_score_r, file=fp)


#########
##################################################
# Important interface function
##################################################
#given multiple sequence alignment
#this funtion will
#return effective counts, target frequencies, and PSSM
def get_neffs_Qs_and_pssm(msa_seqs):

    neffs = []
    feffs = []
    Qs = []
    pssm = []
    for i in range(len(msa_seqs[0])):
        neffs.append({})
        feffs.append({})
        Qs.append({})
        pssm.append({})

    #calculate neffs
    for i in range(len(msa_seqs[0])):
        for a in aa_list2:
            neffs[i][a] = calculate_neff(msa_seqs, i, a)
            #print 'neff',i,a,':',neffs[i][a]

    for i in range(len(msa_seqs[0])):
        #modified not to include gap numbers
        #neff_sum = sum( neffs[i].itervalues() )
        neff_sum = 0
        for k, v in neffs[i].items():
            if k == '-':
                continue
            neff_sum += v

        #print i, neff_sum;
        for a in aa_list:
            if neff_sum:
                feffs[i][a] = neffs[i][a] / neff_sum
            else:
                #feffs[i][a] = 0.0
                #skipping all gap column
                continue

    Nc = calculate_mean_number_of_symbols2(neffs)
    for i in range(len(msa_seqs[0])):
        #skipping all gap column
        if not feffs[i]:
            continue
        for a in aa_list:
            Qs[i][a] = calculate_target_frequency(
                msa_seqs, feffs, i, a, Nc - 1)

    for i in range(len(msa_seqs[0])):
        #skipping all gap column
        if not Qs[i]:
            continue

        for a in aa_list:
            if Qs[i][a] == 0.0:
                #suppose not to print out
                #since the speciall Qs[i][a] == 0.0 should be delt by
                #skipping thing....
                print('WARNING: Q is zero in def get_neffs_Qs_and_pssm!', i, a, file=sys.stderr)
                continue
            try:
                pssm[i][a] = math.log(Qs[i][a] / p_rbnsn[a])
            except:
                print("Error: in calculating pssm", i, a, Qs[i][a], p_rbnsn[a])
                sys.exit()

    return neffs, Qs, pssm

#########################
# save python numerical format as binary format
#########################
#Not implemented correctly!!!!
#########################
#Do not use!!


def save_neffs_Qs_pssm_as_binary_file(neffs, Qs, pssms, fn, overwrite=0):
    if not overwrite:
        if os.path.exists(fn):
            print("Error! %s already exists cannot overwrite!" % fn, file=sys.stderr)
            return

    if len(neffs) == len(Qs) == len(pssms):
        pass
    else:
        print("Error! Neffs, Qs, and PSSMs does not have same dimension!", file=sys.stderr)
        return

    fout = open(fn, 'w')
    lenaa = len(aa_list2)

    ##################
    #Save neffs
    temp = array.array('f', [0.0] * lenaa * len(neffs))
    for i, neff in enumerate(neffs):
        for j, aa in enumerate(aa_list2):  # contains gap (-)
            temp[i * lenaa + j] = neff[aa]
    temp.tofile(fout)

    ###################
    #Save
    temp = array.array('f', [0.0] * lenaa * len(neffs))
    for i, Q in enumerate(Qs):
        for j, aa in enumerate(aa_list2):  # contains gap (-)
            temp[i * lenaa + j] = Q[aa]
    temp.tofile(fout)

    temp = array.array('f', [0.0] * lenaa * len(neffs))
    for i, pssm in enumerate(pssms):
        for j, aa in enumerate(aa_list2):  # contains gap (-)
            temp[i * lenaa + j] = pssm[aa]
    temp.tofile(fout)


def save_neffs_as_binary_file(neffs, fn, overwrite=0):
    if not overwrite:
        if os.path.exists(fn):
            print("Error! %s already exists cannot overwrite!" % fn, file=sys.stderr)
            return

    if len(neffs):
        pass
    else:
        print("Error! Neffs, Qs, and PSSMs does not have same dimension!", file=sys.stderr)
        return

    fout = open(fn, 'w')
    lenaa = len(aa_list2)

    ##################
    #Save neffs
    temp = array.array('f', [0.0] * lenaa * len(neffs))
    for i, neff in enumerate(neffs):
        for j, aa in enumerate(aa_list2):  # contains gap (-)
            temp[i * lenaa + j] = neff[aa]
    temp.tofile(fout)


#########################
#########################
#  MAIN
#########################
#########################

if (__name__ == '__main__'):
    import os
    import pickle
    import glob
    #psyco.full()

    data_dir = sys.argv[1]
    if len(sys.argv) <= 2:
        aln_list = glob.glob(data_dir + '/*.aln')
    else:
        aln_list = sys.argv[2:]

    Qs_list = {}
    neffs_list = {}
    for fn in aln_list:
        id = fn.split('/')[-1][:4]
        print(fn)
        msa_fp = open(fn)

        msa_seqs = parse_msa(msa_fp)
        print(msa_seqs[0])

        neffs = []
        feffs = []
        Qs = []
        pssm = []
        for i in range(len(msa_seqs[0])):
            neffs.append({})
            feffs.append({})
            Qs.append({})
            pssm.append({})

        #calculate neffs
        for i in range(len(msa_seqs[0])):
            for a in aa_list2:
                neffs[i][a] = calculate_neff(msa_seqs, i, a)
                #print 'neff',i,a,':',neffs[i][a]

        for i in range(len(msa_seqs[0])):
            #modified not to include gap numbers
            #neff_sum = sum( neffs[i].itervalues() )
            neff_sum = 0
            for k, v in neffs[i].items():
                if k == '-':
                    continue
                neff_sum += v

            #print i, neff_sum;
            for a in aa_list:
                if neff_sum:
                    feffs[i][a] = neffs[i][a] / neff_sum
                else:
                    feffs[i][a] = 0.0

        Nc = calculate_mean_number_of_symbols2(neffs)
        for i in range(len(msa_seqs[0])):
            for a in aa_list:
                Qs[i][a] = calculate_target_frequency(
                    msa_seqs, feffs, i, a, Nc - 1)

        for i in range(len(msa_seqs[0])):
            for a in aa_list:
                if Qs[i][a] == 0.0:
                    print('WARNING: Q is zero!', i, a)
                    continue
                try:
                    pssm[i][a] = math.log(Qs[i][a] / p_rbnsn[a])  # calculate_pssm( Qs[i][a], p_rbnsn[a] )
                except:
                    print("Error: in calculating pssm", i, a, Qs[
                        i][a], p_rbnsn[a])
                    sys.exit()
            #print pssm[i]

        save_fn = 'junk.dump'
        fp = open(save_fn, 'w')
        for i in range(len(msa_seqs[0])):
            print(i + 1, ':', end=' ', file=fp)
            for a in aa_list2:
                print(neffs[i][a], end=' ', file=fp)
            print(file=fp)

        for i in range(len(msa_seqs[0])):
            print(i + 1, ':', end=' ', file=fp)
            for a in aa_list:
                try:
                    print(pssm[i][a], end=' ', file=fp)
                except KeyError:
                    pass
            print(file=fp)

        #dump_neffs = cPickle.dump( neffs,fp,-1 )
        #dump_pssm = cPickle.dump( pssm,fp,-1 )
        Qs_list[id] = pssm
        neffs_list[id] = neffs

    '''
    neffs = neffs_list['0228']
    pssm = Qs_list['0228']
    fp = open( 'junk.dump','w' )
    for i in xrange( len(msa_seqs[0]) ) :
            print >>fp, i+1,':',
            for a in aa_list2 :
                    print >>fp, neffs[i][a],
            print >>fp

    for i in xrange( len(msa_seqs[0]) ) :
            print >>fp, i+1,':',
            for a in aa_list :
                    print >>fp, pssm[i][a],
            print >>fp

    #cPickle.dump( neffs_list, fp )
    #cPickle.dump( Qs_list, fp )
    '''
