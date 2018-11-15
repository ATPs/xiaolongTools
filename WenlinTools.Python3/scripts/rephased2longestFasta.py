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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa_2list(fa):
    alist = []
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        seq = seq.replace('N', '-')
        alist.append(seq)
    return alist, len(seq)


def read_letter_info(fn):
    #3614_mito	107	T	C	63	3614_mito_1
    rdict = {}
    for line in cmn.file2lines(fn):
        if line.strip() == '':
            continue

        items = line.strip().split()
        scaf, index, char1, char2 = items[:4]
        phase = items[-1]
        #TODO:currently, assume only mito is present
        rdict[int(index) - 1] = (char1, char2, phase)
    return rdict

def find_longest_contig(adict):
    indexes = list(adict.keys())
    indexes.sort()
    i = 0
    j = 0
    current_phase = None
    allcontigs = []
    for index in indexes:
        char1, char2, phase = adict[index]
        if current_phase == None:
            current_phase = phase
            i = index
            j = index

        else:
            #has current_phase
            if current_phase == phase:
                j = index
            else: #phase changed
                allcontigs.append((i, j))
                i = index
                j = index
                current_phase = phase

    allcontigs.append((i, j))

    i, j = max(allcontigs, key=lambda x: x[1] - x[0])
    return i, j


def read_correct_info(fn):
    #3614_mito       5       C
    rdict = {}
    for line in cmn.file2lines(fn):
        items = line.strip().split('\t')
        #TODO: asuming only working on one scaffold (mito)
        index, char, info = items[-3:]
        rdict[int(index) - 1] = char
    return rdict


if __name__=='__main__':
    #options=parse_options()
    try:
        fletter, fa = sys.argv[1:]
    except:
        print("Usage: *.py phased_15103D06_snp_phased_rephased.letters phased_15103D06_snp_phased.fa", file=sys.stderr)
        sys.exit()


    #change 1-base index into 0-base
    letter_dict = read_letter_info(fletter)

    outlabel = '.'.join(fletter.split('.')[:-1]).replace('_rephased', '')

    if len(letter_dict) == 0:
        print('no phasing letters found!')
        dn = outlabel + '_rephased_alltodo.fa'
        cmn.run('cp %s %s' % (fa, dn))
        sys.exit()


    longest_contig_i, longest_contig_j = find_longest_contig(letter_dict)


    sp = cmn.lastName(outlabel).split('_')[1]

    fcorrect = outlabel + '_corrected_positions.txt'

    corrected_dict = read_correct_info(fcorrect)

    seqList, length = read_fa_2list(fa)

    print(outlabel, longest_contig_j, longest_contig_i, length)
    if longest_contig_j - longest_contig_i < 0.5 * length:
        print('less than 50%% are phased for %s, not very useful' % fa)

    newSeq = []
    todoSeq = []
    for i in range(length):
        if i in letter_dict:
            #because the longest contig only cover the phased letters
            if longest_contig_i <= i <= longest_contig_j:
                char1, char2, phase = letter_dict[i]
            else:
                #not in longest contig but has phase
                char1, char2 = 'N', 'N'
        elif i in corrected_dict:
            char = corrected_dict[i]
            char1 = char
            char2 = char

        else:
            char1 = seqList[0][i]
            char2 = seqList[1][i]
            if char1 != char2:
                print('warning: different chars in %s' % i)

        newSeq.append((char1, char2))


        #making todoDict
        #todoseq would not mask the regions
        if i in corrected_dict:
            char = corrected_dict[i]
            char1 = char
            char2 = char
        elif i in letter_dict:
            char1, char2, phase = letter_dict[i]
        else:
            char1 = seqList[0][i]
            char2 = seqList[1][i]

        todoSeq.append((char1, char2))



    dn = outlabel + '_rephased.fa'
    fasta = '>%s_phase1\n%s\n' % (sp, ''.join([each[0] for each in newSeq]))
    fasta += '>%s_phase2\n%s\n' % (sp, ''.join([each[1] for each in newSeq]))

    cmn.write_file(fasta, dn)

    dn = outlabel + '_rephased_alltodo.fa'
    fasta = '>%s_phase1\n%s\n' % (sp, ''.join([each[0] for each in todoSeq]))
    fasta += '>%s_phase2\n%s\n' % (sp, ''.join([each[1] for each in todoSeq]))

    cmn.write_file(fasta, dn)


