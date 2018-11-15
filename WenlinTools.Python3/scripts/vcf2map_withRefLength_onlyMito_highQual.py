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
from random import randint


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_length_info(fn):
    adict = {}
    for line in cmn.file2lines(fn):
        scaf, length = line.strip().split()
        adict[scaf] = int(length)

    #assume the last scaffold is mito
    if 'mito' not in scaf:
        scaf = ''
    return adict, scaf


def fill_in_missing(scaf, smallI, largeI):
    #fill till the index is the same (the same one is not included)

    global mitoScaf, fp_mito, fp_genome
    if scaf == mitoScaf:
        fp = fp_mito
        filled = 'N\n'
    else:
        filled = 'N\tN\n'
        fp = fp_genome

    while(smallI != largeI):
        fp.write(filled)
        smallI += 1


def append_last_chars(scaf, smallI, length):
    global mitoScaf, fp_mito, fp_genome
    if scaf == mitoScaf:
        fp = fp_mito
        filled = 'N\n'
    else:
        fp = fp_genome
        filled = 'N\tN\n'

    if smallI != length:
        #if the length is not meet, fill to the length
        for _ in range(length - smallI):
            fp.write(filled)
            #taken_lines.append('N\tN')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def tell_base_Quality(Ncount, qual_label):
    #LowQual

    Ncount = list(map(int, Ncount))

    #1. firstly check if any of them has 80% more
    Ntotal = sum(Ncount)
    Ncut = 0.8 * Ntotal
    isDominated = any([each >= Ncut for each in Ncount])
    if isDominated:
        return True #good base
    else:
        #reach here if there are confliction
        if Ntotal < 10 or qual_label == 'LowQual':
            return False
        else:
            return True


def vcf2map_mito_highQual(line):
    #high quality means:
    #1. when there are base that is less than 80% of the total base,
    #   if the coverage is lower than 10 or has lowQual label, mask it
    #TODO: mask all the lowQual regions?
    items = line.strip().split()
    scaf, index = items[:2]
    base, mut = items[3:5]

    qual_label = items[6]

    info = items[-1]
    #info_labels = items[-2].split(':')

    subitems = info.split(':')
    GT = subitems[0]

    #if 'LowQual' in line:
        #not sufficient evident to call anything
        #keep original ones
    #    fasta[0].append(base)
    #    fasta[1].append(base)
    #    continue

    char1 = ''
    if GT == './.':#no coverage
        if items[7] == '.':
            char1 = '-'
        else:
            char1 = 'N'
        #print '%s\t%s' % (char1, char2)
        #return [char1, char2]
        return char1

    GT_counts = subitems[1].split(',')

    if GT == '0/0': #no SNP call
        char1 = base
    elif GT == '1/1': #totally different
        char1 = mut

    elif GT == '0/1': #called SNP
        i, j = list(map(int, GT_counts))
        isGoodQual = tell_base_Quality([i, j], qual_label)
        if isGoodQual:
            if i > j:
                char1 = base
            else:
                # when i = j, will make it mut
                char1 = mut
        else:
            char1 = 'N'

    elif GT == '1/2':
        o, i, j = list(map(int, GT_counts))
        if o != 0:
            print('warning: original frequency is not 0 for:', file=sys.stderr)
            print(line, file=sys.stderr)
        base, mut = mut.split(',')

        isGoodQual = tell_base_Quality([o, i, j], qual_label)
        if isGoodQual:
            if i > j:
                char1 = base
            else:
                # when i = j, will make it mut
                char1 = mut
        else:
            char1 = 'N'
    else:
        print('unrecognized line: %s' % line, file=sys.stderr)

    #return [char1, char2]
    #return '%s\t%s' % (char1, char2)
    return char1



def vcf2map(line):
    items = line.strip().split()
    scaf, index = items[:2]
    base, mut = items[3:5]

    info = items[-1]
    subitems = info.split(':')
    GT = subitems[0]
    #if 'LowQual' in line:
        #not sufficient evident to call anything
        #keep original ones
    #    fasta[0].append(base)
    #    fasta[1].append(base)
    #    continue

    char1, char2 = '', ''
    if GT == './.':#no coverage
        if items[7] == '.':
            char1 = '-'
            char2 = '-'
        else:
            char1 = 'N'
            char2 = 'N'
        #print '%s\t%s' % (char1, char2)
        #return [char1, char2]
        return '%s\t%s' % (char1, char2)

    GT_counts = subitems[1].split(',')

    if GT == '0/0': #no SNP call
        char1 = base
        char2 = base
    elif GT == '1/1': #totally different
        char1 = mut
        char2 = mut

    elif GT == '0/1': #called SNP
        i, j = list(map(int, GT_counts))
        randomN = randint(0, 1)
        char1 = base
        char2 = mut
        if randomN:
            char1, char2 = char2, char1

    elif GT == '1/2':
        o, i, j = list(map(int, GT_counts))
        if o != 0:
            print('warning: original frequency is not 0 for:', file=sys.stderr)
            print(line, file=sys.stderr)
        char1, char2 = mut.split(',')
        randomN = randint(0, 1)
        if randomN:
            char1, char2 = char2, char1

    else:
        print('unrecognized line: %s' % line, file=sys.stderr)

    #return [char1, char2]
    return '%s\t%s' % (char1, char2)



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fdict = sys.argv[1:3]
    except:
        print("Usage: *.py *.vcf length.txt", file=sys.stderr)
        sys.exit()

    dnlabel = cmn.lastName(fn).replace('.vcf', '')

    scaf_length, mitoScaf = read_length_info(fdict)

    #used to record the outputted index
    scaf_counters = {}

    #used to make sure the last scaffold output to the end
    lastIndex = -1
    lastScaf = ''

    fgenome = '%s.map' % dnlabel
    fp_genome = open(fgenome, 'w')

    if mitoScaf != '':
        fmito = '%s_highQualMITO.map' % dnlabel
        print('mitogenome detected, would be output into %s' % fmito)
        fp_mito = open(fmito, 'w')



    with open(fn) as fp:
        for line in fp:
            #if not line.startswith('scaffold'):
            if line.strip()[0] == '#':
                continue

            items = line.strip().split()
            scaf = items[0]

            #only take care of mitogenome
            if scaf != mitoScaf:
                continue

            try:
                expect_index = scaf_counters[scaf]
            except:
                expect_index = 1
                #a new scaffold, check if it output to the end
                try:
                    append_last_chars(lastScaf, lastIndex, scaf_length[lastScaf])
                except:
                    if lastScaf != '':
                        print('Error!: scaffold %s not found!' % lastScaf, file=sys.stderr)

            index = int(items[1])
            #make the expect index run the the current
            fill_in_missing(scaf, expect_index, index)

            #parse current
            if scaf == mitoScaf:
                mapline = vcf2map_mito_highQual(line) + '\n'
                fp_mito.write(mapline)
            else:
                mapline = vcf2map(line) + '\n'
                fp_genome.write(mapline)

            scaf_counters[scaf] = index + 1
            lastScaf = scaf
            lastIndex = index


    #append for the last scaffold in the list
    append_last_chars(lastScaf, lastIndex, scaf_length[lastScaf])

    fp_genome.close()
    if mitoScaf != '':
        fp_mito.close()
