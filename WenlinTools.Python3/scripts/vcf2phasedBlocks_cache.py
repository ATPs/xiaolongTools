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
    lengthDict = {}
    lastScaf = None
    lastPhase = None
    with open(fn) as fp:
        for line in fp:
            line = line.strip()
            if line[0] == '#':
                #TODO:get the length dict from the header
                continue

            count += 1
            items = line.strip().split('\t')
            scaf = items[0]
            index = int(items[1])
            base, mut = items[3:5]

            if scaf != lastScaf:
                #output the last scaf
                if lastScaf != None:
                    phased_blocks.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lastPhase, lastScaf, lastPosition[1], right[1], lastPosition[2], right[2]))

                lastScaf = scaf
                #reset the poiters for every scaffold
                lastPosition = [scaf, 1, count, None]
                left = [scaf, 1, count]
                right = [scaf, 1, count]

            infoHeader, info = items[-2:]
            infoHeader_list = infoHeader.split(':')

            info_list = info.split(':')
            GT = info_list[0]
            if GT == './.':#no coverage
                if items[7] == '.':
                    char1 = '-'
                    char2 = '-'
                else:
                    char1 = 'N'
                    char2 = 'N'
                right[1:3] = [index, count]

            elif GT == '0/0': #no SNP call
                char1 = base
                char2 = base
                right[1:3] = [index, count]

            elif GT == '1/1': #totally different
                char1 = mut
                char2 = mut
                right[1:3] = [index, count]

            elif GT == '0/1': #called SNP
                if 'HP' in infoHeader_list:
                    #phased
                    phase1, phase2 = info_list[infoHeader_list.index('HP')].split(',')

                    lastPhase = lastPosition[-1]
                    currentPhase = phase1.split('-')[0]
                    if lastPhase == None:
                        #if only begin, just initiate
                        lastPosition[-1] = currentPhase
                        left[1] = index + 1
                        left[2] = count + 1
                        lastPhase = currentPhase
                    else:
                        if lastPhase == currentPhase:
                            #just continue the phase
                            left[2] = count + 1
                            left[1] = index + 1
                            right[1:3] = (index, count)
                        else:
                            #found a new phase
                            right[1:3] = (index, count)

                            ###1. output
                            phased_blocks.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lastPhase, lastScaf, lastPosition[1], right[1], lastPosition[2], right[2]))
                            #phased_blocks.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lastPhase, scaf, index, lastPosition[1], right[1], lastPhase))

                            ###2, reset
                            lastPosition = list(left)
                            lastPosition.append(currentPhase)

                            ###move left
                            left[1] = index + 1
                            left[2] = count + 1

                    if '-1' in phase1:
                        char1 = base
                        char2 = mut
                    elif '-2' in phase1:
                        char1 = mut
                        char2 = base
                    else:
                        print('Error!!! can not understand this line: %s' % line)
                        sys.exit()
                else:
                    #if this position is not phasable, make it to the dominated letter
                    copy_counts = list(map(int, info_list[infoHeader_list.index('AD')].split(',')))
                    char1, char2 = base, mut
                    if copy_counts[0] > copy_counts[1]:
                        char2 = char1
                    else:
                        char1 = char2

            elif GT == '1/2':
                char1, char2 = mut.split(',')
                if 'HP' in infoHeader_list:
                    phase1, phase2 = info_list[infoHeader_list.index('HP')].split(',')

                    lastPhase = lastPosition[-1]
                    currentPhase = phase1.split('-')[0]
                    if lastPhase == None:
                        #if only begin, just initiate
                        lastPosition[-1] = currentPhase
                        left[2] = count + 1
                        left[1] = index + 1
                        lastPhase = currentPhase
                    else:
                        if lastPhase == currentPhase:
                            #just continue the phase
                            left[2] = count + 1
                            left[1] = index + 1
                            right[1:3] = [index, count]
                        else:
                            #found a new phase
                            right[1:3] = [index, count]

                            ###1. output
                            phased_blocks.append('%s\t%s\t%s\t%s\t%s\n' % (lastPhase, scaf, index, lastPosition[1], right[1], lastPhase))

                            ###2, reset
                            lastPosition = list(left)
                            lastPosition.append(currentPhase)

                            ###move left
                            left[1] = index + 1
                            left[2] = count + 1

                    if '-1' in phase1:
                        pass
                    if '-2' in phase1:
                        char1, char2 = char2, char1
                    else:
                        print('Error!!! can not understand this line: %s' % line)
                        sys.exit()
                else:
                    #if this position is not phasable, make it to the dominated letter
                    copy_counts = list(map(int, info_list[infoHeader_list.index('AD')].split(',')))[1:]
                    if copy_counts[0] > copy_counts[1]:
                        char2 = char1
                    else:
                        char1 = char2

            else:
                print('unrecognized line: %s' % line, file=sys.stderr)
                sys.exit()


            seq1.append(char1)
            seq2.append(char2)

    #output the last
    phased_blocks.append('%s\t%s\t%s\t%s\t%s\t%s\n' % (lastPhase, lastScaf, lastPosition[1], right[1], lastPosition[2], right[2]))

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


    dn = dnlabel + 'phased.blocks'
    cmn.write_file(''.join(phased_blocks), dn)
