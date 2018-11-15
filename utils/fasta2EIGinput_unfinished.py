
import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os
from collections import Counter

def fasta2EIGinput(filename,gapcut = 0.8, Ncut=4):
    '''
    given the filename of a fasta file, generate files required to run PCA analysis
    filename, input fasta file
    gapcut, keep sites with base ratio larger than gapcut, default 0.8
    Ncut, count SNP only if the SNP exist more than this number. Default 4
    '''
    filename_base = os.path.basename(filename)
    filename_dir = os.path.dirname(filename)
    

if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    gapcut = 0.8 #the positions has to have more than this number of chars
    #enrich more informative positions
    Ncut = 4 #each nucleatide appear should have larger count than this

    #fnlabel = cmn.lastName(fn) + '_loose'
    fnlabel = cmn.lastName(fn)

    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].split('_')[0]
            else:
                seq = line.strip()
                try:
                    adict[name].append(seq)
                except KeyError:
                    adict[name] = [seq]
    length = len(seq)
    #snp file like 'No. chr. Mdist pos char1 char2'
    #here I put everything as one chromosome to treat the program
    #and chromsome number start at 100 because lower number has specical assignment in program

    #the program only accept biallelic positions

    #ind file like 'sample U Ignore'

    #geno file like '012' meaning how many of char1 appeared

    names = list(adict.keys())

    #make ind file
    #TODO: make pop assignment
    ind_lines = ['%s\tU\tcase1\n' % (name) for name in names]
    dn = fnlabel + '.ind'
    cmn.write_file(''.join(ind_lines), dn)

    snp_lines = []
    geno_lines = []

    missing_data = set('N -'.split())
    #header = ['fake', 'position'] + names
    #header = '\t'.join(header)


    count = 1
    for index in range(length):
        char_dict = {name:[seq[index] for seq in adict[name]] for name in adict}
        chars = sum(list(char_dict.values()), [])
        check_chars = [char for char in chars if (char != '-') and (char != 'N')]

        #filter gaps
        if len(check_chars) < (gapcut * len(chars)):
            continue

        #filtering informative position
        check_chars_set = set(check_chars)
        N = len(check_chars_set)
        count_dict = Counter(check_chars)

        goodChars = [char for char in count_dict
                if count_dict[char] >= Ncut]
        check_chars_set = set(goodChars)
        N = len(check_chars_set)
        count_dict = Counter(check_chars)


        if N == 2: #only consider bialletic positions

            ref_char = max(list(count_dict.keys()), key=lambda x: count_dict[x])
            other_char = [char for char in check_chars_set if char != ref_char][0]
            count += 1

            #snp file like 'No. chr. Mdist pos char1 char2'
            snp_line = ['snp%s' % count, 10, '0.0', index, ref_char, other_char]
            snp_lines.append(' '.join(map(str, snp_line)))

            geno_line = []
            for i, name in enumerate(names):
                spchars = char_dict[name]
                if any([char in missing_data for char in spchars]):
                    geno_line.append(9)
                    continue

                code = spchars.count(ref_char)
                geno_line.append(code)
            geno_lines.append(''.join(map(str, geno_line)))

    dn = fnlabel + '.geno'
    geno_lines.append('')
    cmn.write_lines(geno_lines, dn)

    dn = fnlabel + '.snp'
    snp_lines.append('')
    cmn.write_lines(snp_lines, dn)


