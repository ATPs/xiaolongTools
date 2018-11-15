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
import pysam
from collections import Counter

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def find_insertion_in_ref(aligned, seq):
    index = 1
    startIndel = False
    adict = {}
    letters = []
    leftI = None
    rightI = None
    lastJ = aligned[0][1]
    while(index < len(aligned) - 1):
        i, j = aligned[index]

        if j == None and i != None:
            if startIndel:
                #extending an indel
                letters.append(seq[i])
            elif lastJ != None:
                #found indel
                leftI = lastJ
                rightI = lastJ + 1
                startIndel = True
                letters.append(seq[i])
        elif j != None and startIndel:
            startIndel = False
            adict[(leftI, rightI)] = ''.join(letters)
            letters = []
        lastJ = j
        index += 1

    return adict

def read_baits():
    fns = cmn.cmd2lines('ls baits/bait*.fa')
    seqDict = {}
    for fn in fns:
        name, seq = cmn.file2lines(fn)
        seqDict[name[1:]] = list(seq)
    return seqDict

def update_baits(bait_dict):
    adict = {}
    for i, name in enumerate(bait_dict):
        fnlabel = 'bait%s' % i
        dn = 'baits/%s.fa' % fnlabel
        seq = bait_dict[name]
        fasta = '>%s\n%s\n' % (name, ''.join(seq))
        cmn.write_file(fasta, dn)
        cmd = 'module add bwa; bwa index %s -p %s' % (dn, fnlabel)
        cmn.run(cmd)
        adict[name] = dn
    return adict

def group_fq(fqlist):
    adict = {}
    for fn in fqlist:
        fnlabel = cmn.lastName(fn)
        sp = fnlabel.split('_')[0]
        if sp not in adict:
            adict[sp] = [None, None, None]
        if '_R1' in fnlabel:
            adict[sp][0] = fn
        elif '_R2' in fnlabel:
            adict[sp][1] = fn
        elif '_singleton' in fnlabel:
            adict[sp][2] = fn
        else:
            print('Error! Can not recognize file %s' % fn)
    return adict

if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py samlist", file=sys.stderr)
        sys.exit()

    fnlist = cmn.file2lines(fn)
    pCoverage = {}

    indel_dict = {}

    for fn in fnlist:
        if not cmn.filexist(fn):
            continue
        print('checking %s ...' % fn)
        for record in pysam.AlignmentFile(fn):
            if record.is_unmapped or record.is_secondary:
                continue
            #i is read, j is ref
            aligned = record.aligned_pairs
            read_seq = record.query_sequence
            for i, j in aligned:
                try:
                    pCoverage[j] += 1
                except KeyError:
                    pCoverage[j] = 1


            #if j == None and i != None:
            #find insertions in reference
            indels = find_insertion_in_ref(aligned, read_seq)
            #print record.query_name, indels, aligned
            for key in indels:
                try:
                    indel_dict[key].append(indels[key])
                except KeyError:
                    indel_dict[key] = [indels[key]]

    bait_dict = read_baits()
    bait_names = list(bait_dict.keys())

    cmn.run('rm bait_insertion 2> /dev/null')
    hasInsertion = False
    for key in indel_dict:
        print(key, len(indel_dict[key]))
        leftI, rightI = key
        cov = (pCoverage[leftI] + pCoverage[rightI]) / 2.0
        indel_depth = len(indel_dict[key])
        insertion_info = [key , indel_depth, cov]
        if indel_depth > 0.5 * cov:
            hasInsertion = True
            count_dict = Counter(indel_dict[key])
            maxChar = max(count_dict, key=lambda x: count_dict[x])
            insertion_info.append(maxChar)
            insertion_info = '\t'.join(map(str, insertion_info))
            cmn.append_file(insertion_info + '\n', 'bait_insertion')

            print('insert between ', insertion_info)
            for name in bait_names:
                bait_dict[name][leftI] += maxChar

    #just undergo one round of adding gap
    if not hasInsertion:
        print('No need to re-run bwa because no insertion in query')
    else:
        N = cmn.cpu_check()
        print('re-run bwa due to insertion')
        frefs = update_baits(bait_dict)

        fq_groups = group_fq(cmn.file2lines('fqlist'))

        bwa_cmds = []
        for reflabel in frefs:
            fref = frefs[reflabel]
            fnlabel = cmn.lastName(fref).replace('.fa', '')
            for sp in fq_groups:
                R1, R2, single = fq_groups[sp]
                cmd = 'bwa mem -B 2 -t %s -M %s %s %s | grep "%s" > %s_paired_%s_mapped.sam ' % (N, fnlabel, R1, R2, reflabel, sp, fnlabel)
                bwa_cmds.append(cmd)
                cmd = 'bwa mem -B 2 -t %s -M %s %s | grep "%s" > %s_single_%s_mapped.sam ' % (N, fnlabel, single, reflabel, sp, fnlabel)
                bwa_cmds.append(cmd)
                bwa_cmds.append('\nwait\n')

        dn = 'bwa1_rerun.cmds'
        cmn.write_lines(bwa_cmds, dn)
        cmd = 'bash %s > insertionRerun.log' % dn
        cmn.run(cmd)





