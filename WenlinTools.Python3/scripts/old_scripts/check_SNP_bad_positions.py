#!/usr/bin/env python

#function:identify bad positions in SNP calls and then report the read coverage
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home2/wli/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
from collections import Counter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def load_alignment_2dict(fn):
    fastas = cmn.txt_read(fn).split('>')[1:]

    adict = {}
    for fasta in fastas:
        lines = fasta.split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        adict[defline] = seq
    return adict


def check_bad_chars(chars, check_cutoff):
    bad_chars = set([])
    count_dict = Counter(chars)
    for char in count_dict:
        count = count_dict[char]
        if count <= check_cutoff:
            bad_chars.add(char)
    if '-' in bad_chars:
        bad_chars.remove('-')
    return bad_chars


def make_char_stat(chars):
    #bases = 'A T G C -'.split()
    adict = Counter(chars)
    values = [(key, adict[key]) for key in adict]
    sorted_values = sorted(values, cmp=lambda x, y: cmp(x[1], y[1]), reverse=True)

    stat = '|'.join(['%s:%s' % value for value in sorted_values])
    return stat

def process_vcf(fn):
    adict = {}
    with open(fn) as fp:
        for line in fp:
            #scaffold1_cov14552_reverse  52  .   C   T   207.74  .   AC=1;AF=0.500;AN=2;BaseQRankSum=0.348;DP=10;Dels=0.00;FS=3.979;HaplotypeScore=0.000
            #0;MLEAC=1;MLEAF=0.500;MQ=57.34;MQ0=0;MQRankSum=1.393;QD=20.77;ReadPosRankSum=1.393;SOR=1.198 GT:AD:DP:GQ:PL  0/1:1,9:10:6:235,0,6

            #scaffold1_cov14552_reverse	20	.	A	.	40.23	.	AN=2;DP=4;MQ=60.00;MQ0=0	GT:DP	0/0:4

            if line.startswith('#'):
                continue
            items = line.strip().split('\t')
            index = int(items[1])
            called = items[4]
            if called == '.':
                called = items[3]
            ref = items[3]
            info = items[9]
            key = '%s\t%s' % (index, called)

            if info == './.':#gap position
                stat = '%s\t0\t0\t0' % ref
            else:
                freq = info.split(':')[1]

                try:
                    i, j = list(map(int, freq.split(',')))
                    total = i+j
                    stat = '%s\t%s\t%s\t%.2f' % (ref, j, total, (j/float(total)))
                except:
                    if len(freq.split(',')) == 1:
                        stat = '%s\t%s\t%s\t1.0' % (ref, freq, freq)
                    else:
                        stat = 'specical_frequency'
            adict[key] = stat
    return adict

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #run this in such directory as
    #/project/biophysics/Nick_lab/wli/sequencing/3318_mito/6_build_tree/2_verify_alignment

    f_aln = '../1_process_alignment/all_genomes.fa'

    fasta_dict = load_alignment_2dict(f_aln)

    #find the bad positions
    check_cutoff = 6 #if less than 6 appear in the count, need to check this position

    #logic: for positions, check for the minor characters
    #logic: if the characters <= check_cutoff, need to go to reads to find its coverage
    #logic: I don't check gap: gap means 0 coverage, just report them
    #alg: firstly find such position for each sequence, then go to their corresponding vcf files

    #TODO: check coverage
    keys = list(fasta_dict.keys())
    seqs = [fasta_dict[key] for key in keys]
    #pos_dict:positions to be checked
    pos_dict = {key:[] for key in keys}
    #gap_dict: positions of gaps
    gap_dict = {key:[] for key in keys}
    #stat_dict: record the frequency of a position
    stat_dict = {}

    for position in range(len(seqs[0])):
        #Note: this position starts with 0
        #so need the good position to +1
        chars = [seq[position] for seq in seqs]
        if len(set(chars)) == 1:
            continue

        #bad_chars = check_bad_chars(chars, check_cutoff)
        bad_chars = set(chars)
        stat_dict[position+1] = make_char_stat(chars)
        for index, char in enumerate(chars):
            if char in bad_chars:
                pos_dict[keys[index]].append('%s\t%s' % (position+1, char))
            if char == '-':
                gap_dict[keys[index]].append(position+1)

    header = 'position called_base ref_base countN totalN percentage stat_in_aln'.split()
    for key in keys:
        new = ['\t'.join(header)]
        print('processing %s' % key)
        label = key.split('_')[0]
        positions = pos_dict[key]
        f_vcf = '%s_snp.vcf' % label
        #vcf_dict = process_vcf(f_vcf)
        puStat_dict = process_puStat(f_puStat)
        for position in positions:
            i, char = position.split()
            stat = stat_dict[int(i)]
            try:
                freq = vcf_dict[position]
            except:
                freq = 'Need Manual'
            newline = '%s\t%s\t%s' % (position, freq, stat)
            new.append(newline)
        dn = '%s_SNP_check.info' % label
        cmn.write_lines(new, dn)
