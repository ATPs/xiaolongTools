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
def fill_in_missing(smaller_index, index):
    #fill in the difference between index
    N = index - smaller_index
    return 'N' * N


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py vcf", file=sys.stderr)
        sys.exit()


    length_dict = {}

    current_scaf = ''
    seqDict = {}
    order_list = []
    with open(fn) as fp:
        for line in fp:
            line = line.strip()

            # ##contig=<ID=scaffold1_cov51,length=30279>
            if line.startswith('##contig='):
                scaf = cmn.find_between(line, '<ID=', ',')
                length = int(cmn.find_between(line, ',length=', '>'))
                length_dict[scaf] = length


            if line[0] != '#':
                items = line.strip().split()
                scaf = items[0]
                if scaf != current_scaf:
                    order_list.append(scaf)
                    #start a new scaffold
                    expect_index = 1
                    current_scaf = scaf
                else:
                    expect_index += 1

                index = int(items[1])

                char = []
                #print expect_index, index
                if expect_index != index:
                    char += fill_in_missing(expect_index, index)
                    expect_index = index

                char.append(items[3])

                try:
                    seqDict[scaf] += char
                except KeyError:
                    seqDict[scaf] = char


    dn = fn + '_v2r.fa'
    with open(dn, 'w') as dp:
        for scaf in order_list:
            seq = ''.join(seqDict[scaf])
            missingEnd = length_dict[scaf] - len(seq)
            if missingEnd > 0:
                seq = seq + 'N' * missingEnd
            fasta = '>%s\n%s\n' % (scaf, seq)
            dp.write(fasta)
