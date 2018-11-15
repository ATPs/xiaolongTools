#!/usr6/local/bin/python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import cmn

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def split_pollyA(seq):
    #firstly check A position
    #if pollyA can be linked longer than 10, then split it

    alist = []
    i = 0
    j = 0
    for k, char in enumerate(seq):
        if char != 'A':
            if i != j:
                alist.append((i+1, j+1))
            i = k
            j = k
        else:
            j += 1

    #filtering
    tmplist = []
    cutI = 30
    cutJ = len(seq) - 30
    for i, j in alist:
        if j - i < 3:
            continue
        if i < cutI or j > cutJ:
            continue
        tmplist.append((i, j))

    alist = tmplist

    if len(alist) == 0:
        return None, None
    #extend
    blist = []
    i, j = alist[0]
    for i1, j1 in alist[1:]:
        if (i1 - j) <= 1:
            j = j1
        else:
            blist.append((i, j))
            i = i1
            j = j1

    blist.append((i, j))

    longest = max(blist, key=lambda x: x[1] - x[0])
    return longest



def find_match_region(seq1, seq2, misP, matchN):
    #firstly, try to find identical match
    i, j = 0, 0
    ranges = []
    for k in range(len(seq1)):
        j = k
        if seq1[k] != seq2[k] and i != j:
            ranges.append((i, j))
            i = k + 1
            j = k + 1
    if i != j:
        ranges.append((i, j))

    #try to link nearby
    lefts = [each[0] for each in ranges]
    rights = [each[1] for each in ranges]
    #maxN = 0
    for left in lefts:
        for right in rights:
            if left >= right:
                continue
            length = right - left
            #if length <= maxN:
            #    continue
            diffN = compute_hamming_distance(seq1[left:right], seq2[left:right])
            if diffN <= (misP * length):
                #maxN = max(maxN, length)
                if length > matchN:
                    return True
    return False
    #return maxN



def filter_adaptor(seq, ad):
    seedN = 10 #seed used initiate search
    matchN = 25 #minial threshold to consider as a match
    mismatchP = 0.05 # level of allowed mismatch

    seed_M = []
    for k in range(0, len(ad) - seedN):
        seed = ad[k:k+seedN]
        i = seq.find(seed)
        if i != -1:
            seed_M.append((i, k))

    if len(seed_M) == 0:
        #no adaptor match
        return seq, len(seq), False
    else:
        #check if there is any soft match

        #1. get the potential regions
        regions = set([])
        for i, k in seed_M:
            seedi, seedj = 0, len(ad)
            left = (i - k)
            right = left + len(ad) + 1
            if left < 0:
                seedi -= left
                left = 0
            if right > len(seq):
                seedj -= (right - len(seq))
                right = len(seq)

            regions.add((left, right, seedi, seedj))

        #check if any region matches
        matches_left = set([])
        for left, right, adleft, adright in regions:
            if left in matches_left:
                continue

            adRegion = ad[adleft: adright]
            seqRegion = seq[left: right]
            isMatch = find_match_region(adRegion, seqRegion, mismatchP, matchN)
            if isMatch:
                matches_left.add(left)
        if len(matches_left) == 0:
            return seq, len(seq), False
        left = min(matches_left)
        return seq[:left], left, True

def compute_hamming_distance(a, b):
    #print a
    #print b
    seqlength = len(a)
    #if len(b) < len(a):
    #    seqlength = len(b)
    if len(a) != len(b):
        return None

    check_list = [a[i] != b[i] for i in range(seqlength)]
    return sum(check_list)

if __name__=='__main__':
    #options=parse_options()
    try:
        fR1, fR2, ad1, sample = sys.argv[1:]
    except:
        print("Usage: *.py R1 R2 ad1 sample", file=sys.stderr)
        print("ad1 is the adaptor showing up in R1", file=sys.stderr)
        sys.exit()

    #this is already the reverse direction
    ad2 = 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'

    #dn1 = 'filtered_' + cmn.lastName(fR1)
    #dn2 = 'filtered_' + cmn.lastName(fR2)
    dn1 = '%sflt_R1.fastq' % sample
    dn2 = '%sflt_R2.fastq' % sample


    totalline = 0
    totalN = 0
    goodline = 0
    goodN = 0
    with open(fR1) as fp1, open(fR2) as fp2, open(dn1, 'w') as dp1, open(dn2, 'w') as dp2:
        for count, line1 in enumerate(fp1):
            line2 = fp2.readline()
            if count % 4 == 0:
                defline1 = line1
                defline2 = line2
            elif count % 4 == 1:
                seq1 = line1.strip()
                seq2 = line2.strip()
                totalN +=  len(seq1) + len(seq2)
                totalline += 1
                filtered_seq1, left1, isMatch1 = filter_adaptor(seq1, ad1)
                filtered_seq2, left2, isMatch2 = filter_adaptor(seq2, ad2)
                #print left1, left2
                if filtered_seq1 != '' and filtered_seq2 != '':
                    if isMatch1 and isMatch2:
                        print('Warning: sequence cut before adaptor for %s' % defline1.split()[0])
                    dp1.write(defline1)
                    dp2.write(defline2)
                    dp1.write(filtered_seq1+'\n')
                    dp2.write(filtered_seq2+'\n')
                    goodline += 1
                    goodN += len(filtered_seq1) + len(filtered_seq2)
                    isGood = True
                else:
                    isGood = False

            elif count % 4 == 2:
                if isGood:
                    dp1.write(line1)
                    dp2.write(line2)
            elif count % 4 == 3:
                if isGood:
                    dp1.write(line1[:left1]+'\n')
                    dp2.write(line2[:left2]+'\n')

    stat = '\t'.join(map(str, [dn1, dn2, totalline, goodline, totalN, goodN, float(goodline)/totalline, float(goodN)/totalN]))
    #dn = cmn.lastName(fR1).split('_')[0] + 'filtering.stat'
    dn = '%s_filtering.stat' % sample
    cmn.write_file(stat, dn)
