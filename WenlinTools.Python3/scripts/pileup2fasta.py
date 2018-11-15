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

def parse_pileup_to_chars(fn):
    adict = {}
    with open(fn) as fp:
        for line in fp:
            items = line.strip().split()
            #scaffold9_len908707_cov91       23050   N       10      GGggGggg^6G^6g  JJJJJJJJAJ
            key = '%s-%s' % (items[0], items[1])
            letters = items[4]
            char = get_dominated_char(letters)
            adict[key] = char
    return adict

def get_dominated_char(letters):
    count_dict = {}
    for char in 'ATCG':
        count_dict[char] = letters.count(char)

    totalN = sum(count_dict.values())
    for char in count_dict:
        if count_dict[char] > (0.8 * totalN):
            return char

    return 'N'


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, fhead, outlabel = sys.argv[1:]
    except:
        print("Usage: *.py fpileup fheader outlabel", file=sys.stderr)
        sys.exit()



    print('reading pileup file...')
    char_dict = parse_pileup_to_chars(fn)

    seq = []

    print('making fasta...')
    with open(fhead) as fp:
        for i, line in enumerate(fp):
            scaf, index = line.strip().split()
            key = '%s-%s' % (scaf, index)
            try:
                char = char_dict[key]
            except KeyError:
                char = '-'
            seq.append(char)

    dn = '%s.fa' % outlabel
    fasta = '>%s\n%s\n' % (outlabel, ''.join(seq))
    cmn.write_file(fasta, dn)


