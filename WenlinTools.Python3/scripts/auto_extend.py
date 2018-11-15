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
import os
from collections import Counter


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rdict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        '-': 'N',
        '*': 'N'
        }


def grep_reads(read, f_libs, direction):
    #reverse the read
    #reverse = ''.join([rdict[i] for i in read[::-1]])

    cmds = []
    for fn in f_libs:
        cmd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/grep_reads.py %s %s %s &' % (read, fn, direction)
        cmds.append(cmd)

    cmds.append('\nwait;\n')

    f_job = 'grep_read.job'
    cmn.write_lines(cmds, f_job)

    cmn.run('bash %s ' % f_job)
    #the output dir is grep_out
    dn = 'all_grep_reads.txt'
    cmn.run('cat grep_out/* > %s' % dn)
    fished_reads = cmn.getid(dn)

    return fished_reads


def detect_dominated_reads(seqs, seed_i):
    global fraction
    global minN

    sequence_index = seed_i
    form_dict = {} #record which seq has how many faction
    forms = [] #used for output
    while (True):
        subseqs = [seq[:sequence_index] for seq in seqs
                if len(seq) >= sequence_index]
        if len(subseqs) == 0:
            break
        total = float(len(subseqs))

        adict = Counter(subseqs)
        #print 'seq', seqs
        #print 'sequence_index', sequence_index
        #print 'subseq', subseqs
        #print adict

        maxCount = 0
        maxSeq = ''
        for seq in adict:
            count = adict[seq]
            if count > maxCount:
                maxCount = count
                maxSeq = seq

        forms.append('%.3f\t%s\t%s' % (maxCount/total, len(subseqs), maxSeq))
        form_dict[maxSeq] = [maxCount/total, len(subseqs)] #fraction of reads and number of remaining reads

        if maxCount < (fraction * total):
            break
        if len(subseqs) < minN:
            break

        sequence_index += 1

    dn = 'form_stat.txt'
    cmn.write_lines(forms, dn)
    return form_dict

def get_extended(adict):
    global fraction
    global minN

    #kinda very aggresive take the longest possible read

    longest = ''
    for seq in adict:
        f, N = adict[seq]
        #print f, fraction
        #print N, minN
        if f >= fraction and N >= minN:
            if len(seq) > len(longest):
                longest = seq
        #print longest
    return longest

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        seed_read, direction, outlabel, fqlist = sys.argv[1:5]
    except:
        print("Usage: *.py seed direction outdir fqlist", file=sys.stderr)
        print("Usage: *.py AATTAATTAGAAATGAAATGTTAATCGTTTTTAAAAATAT forward 3318 fqlist [100]", file=sys.stderr)
        sys.exit()

    ###################################
    #parameters
    #1. overlaps for extension
    overlapN = 40 # somewhat agresive
    #2. minimal number of reads (coverage) to trust the extension
    try:
        minN = int(sys.argv[4])
    except:
        minN = 100

    print('minimal coverage is %s' % minN)
    #3. faction cutoff to trust the extension
    fraction = 0.5 #50%
    #4. upper limit for iteration
    upper_iter = 10
    #5. cmd to get the read files
    #cmd = 'ls /project/biophysics/Nick_lab/wli/sequencing/3318_2nd/6.5_repair_ID/*.fastq'
    #cmd = 'ls /project/biophysics/Nick_lab/wli/sequencing/3318_2nd/6.5_repair_ID/3318ac_6KP_R1.fastq'
    #f_libs = [os.path.abspath(i) for i in cmn.cmd2lines(cmd)]
    f_libs = [os.path.abspath(line) for line in cmn.file2lines(fqlist)]

    #####################################

    pdir = '%s_extend' % outlabel
    cmn.mkdir(pdir)
    os.chdir(pdir)

    if direction == 'backward':
        print('reverse the order of sequence (not reverse strand)')
        seed_read = seed_read[::-1]

    read = seed_read
    Iter = 0
    extensions = []
    while(Iter < upper_iter):
        Iter += 1
        print('running iteration %s' % Iter)

        #prepare wdir
        wdir = 'extend_iter%s' % Iter
        cmn.mkdir(wdir)
        os.chdir(wdir)
        cmn.write_file(read, 'seed_seq.txt')

        #grep the reads
        fished_reads = grep_reads(read, f_libs, direction)
        cmn.write_lines(fished_reads, 'fished_reads.txt')

        #make stat of the reads
        stat_dict = detect_dominated_reads(fished_reads, len(read))

        #get the longest set under cutoff
        extended_seq = get_extended(stat_dict)
        cmn.write_file(extended_seq, 'extension_seq.txt')
        os.chdir('..')

        if extended_seq == '':
            print('no extension can be found in iteration %s! exit!' % Iter)
            break
        read = extended_seq[-overlapN:]
        extensions.append(extended_seq)


        if len(extensions) >= 2:
            if extensions[-1] == extensions[-2]:
                #no extension now
                break


    #combine the extensions
    Eread = extensions[0]
    for each in extensions[1:]:
        Eread += each[overlapN:]

    dn = '%s_%s_iter%s.extend' % (outlabel, direction, upper_iter)
    cmn.write_file(Eread, dn)




