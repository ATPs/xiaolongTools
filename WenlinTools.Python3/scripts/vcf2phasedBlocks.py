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
def output_phase_position(phase1, phase2, base, mut):
    if '-1' in phase1:
        char1 = base
        char2 = mut
    elif '-2' in phase1:
        char1 = mut
        char2 = base
    else:
        print('Error!!! can not understand this line: %s' % line)
        sys.exit()
    return char1, char2, phase1.split('-')[0]


def merge_lowQual_position(copy_counts, char1, char2):
    if copy_counts[0] > copy_counts[1]:
        char2 = char1
    else:
        char1 = char2
    return char1, char2


def phase_letter(items):
    isHetero = False #flag to label heterozyocity position
    isPhased = False
    global curPhase

    infoHeader, info = items[-2:]
    infoHeader_list = infoHeader.split(':')
    info_list = info.split(':')
    GT = info_list[0]

    base, mut = items[3:5]

    if GT == './.':#no coverage
        if items[7] == '.':
            char1 = '-'
            char2 = '-'
        else:
            char1 = 'N'
            char2 = 'N'

    elif GT == '0/0': #no SNP call
        char1 = base
        char2 = base

    elif GT == '1/1': #totally different
        char1 = mut
        char2 = mut

    elif GT == '0/1': #called SNP
        isHetero = True
        if 'HP' in infoHeader_list:
            #phased
            isPhased = True
            phase1, phase2 = info_list[infoHeader_list.index('HP')].split(',')
            char1, char2, curPhase = output_phase_position(phase1, phase2, base, mut)

        else:
            #if this position is not phasable, make it to the dominated letter
            copy_counts = list(map(int, info_list[infoHeader_list.index('AD')].split(',')))
            char1, char2 = merge_lowQual_position(copy_counts, base, mut)

    elif GT == '1/2':
        isHetero = True
        char1, char2 = mut.split(',')
        if 'HP' in infoHeader_list:
            isPhased = True
            phase1, phase2 = info_list[infoHeader_list.index('HP')].split(',')
            char1, char2, curPhase = output_phase_position(phase1, phase2, char1, char2)

        else:
            #if this position is not phasable, make it to the dominated letter
            copy_counts = list(map(int, info_list[infoHeader_list.index('AD')].split(',')))[1:]
            char1, char2 = merge_lowQual_position(copy_counts, char1, char2)

    else:
        print('unrecognized line: %s' % line, file=sys.stderr)
        sys.exit()

    return char1, char2, isHetero, isPhased


def initiate_counter(scaf):
    global count
    left = [scaf, 1, count, None]
    lastPosition = list(left)
    right = list(left)
    return lastPosition, left, right

def move_pointer(p, index, count):
    global curPhase
    p[1:4] = [index, count, curPhase]
    return p

def ouput_phased_blocks(left, right, label='NA'):
    global lastPhase, NphaseChar
    scaf, index1, count1, phase1 = left
    scaftmp, index2, count2, phase2 = right

    #line = '%s\t%s\t%s\t%s\t%s\t(phase:%s)\n' % (scaf, index1, index2, count1, count2, lastPhase)
    line = '%s\t%s\t%s\t%s\t%s\t%s\t(phase:%s,%s)\n' % (scaf, index1, index2, count1, count2, count2-count1, lastPhase, NphaseChar)
    NphaseChar = 0
    return line

def parse_scaf_length(line):
    ##contig=<ID=scaffold1_len454574_cov98,length=458273>
    adict = {}
    scaf = cmn.find_between(line, 'ID=', ',length')
    length = int(cmn.find_between(line, ',length=', '>'))
    adict[scaf] = length
    return adict


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    #use three pointers to solve this problem
    #lastPosition is the last phased block
    lastPosition = None #the left most of the phased blocks
    #left would move to the candidates of last phased block
    left = None #the index used to record the left most of a block of identical sequence
    #right is the current phased block
    right = None # the right most pointer
    #each pointer is an object contains (scaf, index_in_scaffold, index_in_fasta(pythonIndex))

    count = -1
    seq1 = []
    seq2 = []
    phased_blocks = []

    lastPhase = None
    curPhase = None #used to label the current phase
    lastScaf = None
    NphaseChar = 0

    phased_letter = []
    len_dict = {}
    expectIndex = 1
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '#':
                #TODO:get the length dict from the header
                ##contig=<ID=scaffold1_len454574_cov98,length=458273>
                if line.startswith('##contig'):
                    len_dict.update(parse_scaf_length(line))
                continue

            count += 1
            items = line.strip().split('\t')
            scaf = items[0]
            index = int(items[1])

            #now assuming every position is either phasable or lowQual, so isHeter flag is not usable here
            char1, char2, isHetero, isPhased = phase_letter(items)
            if isPhased:
                phased_letter.append('%s\t%s\t%s\t%s\t(phase:%s)' % (scaf, index, char1, char2, curPhase))

            #now, we have assumed that every position can be phased
            #not phasable position would be merge to homozyocity

            if scaf != lastScaf:
                if lastScaf != None:
                    #it has to have a very last piece to output
                    phased_blocks.append(ouput_phased_blocks(lastPosition, right, 'init%s' % lastScaf))
                    #check to fill in the last
                    scafLength = len_dict[lastScaf]
                    if scafLength + 1 != expectIndex:
                        gap = scafLength + 1 - expectIndex
                        fillN = 'N' * gap
                        seq1.append(fillN)
                        seq2.append(fillN)

                expectIndex = 1
                cutPhase = None
                lastPhase = None
                lastPosition, left, right = initiate_counter(scaf)
                lastScaf = scaf

            if index != expectIndex:
                gap = index - expectIndex
                count += gap
                fillN = 'N' * gap
                seq1.append(fillN)
                seq2.append(fillN)
                expectIndex = index

            seq1.append(char1)
            seq2.append(char2)
            expectIndex += 1

            #moving the pointers
            right = move_pointer(right, index, count)
            #print 'right', right

            #checking if we have good phase
            if isPhased:
                NphaseChar += 1
                if curPhase != lastPhase:
                    if lastPhase == None:
                        lastPhase = curPhase
                    else:
                        #output
                        phased_blocks.append(ouput_phased_blocks(lastPosition, right))
                        lastPosition = list(left)
                        lastPhase = curPhase
                left = move_pointer(left, index + 1, count + 1)
            #print 'left', left

    scafLength = len_dict[scaf]
    if scafLength + 1 != expectIndex:
        gap = scafLength + 1 - expectIndex
        fillN = 'N' * gap
        seq1.append(fillN)
        seq2.append(fillN)

    #output the last
    phased_blocks.append(ouput_phased_blocks(lastPosition, right, 'lastOne'))

    dnlabel = cmn.lastName(fn).replace('.vcf', '')
    sp = dnlabel.split('_')[1]
    dn = dnlabel + '_phased.fa'
    with open(dn, 'w') as dp:
        dp.write('>%s_ref_or_phase1\n' % sp)
        dp.write(''.join(seq1))
        dp.write('\n')
        dp.write('>%s_called_or_phase2\n' % sp)
        dp.write(''.join(seq2))
        dp.write('\n')


    dn = dnlabel + '_phased.blocks'
    cmn.write_file(''.join(phased_blocks), dn)

    dn = dnlabel + '_phased.letters'
    cmn.write_lines(phased_letter, dn)
