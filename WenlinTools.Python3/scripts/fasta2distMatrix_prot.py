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
dist_records = {}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        defline, seq = each.strip().split()
        seq = list(seq.strip())
        #Ngap = seq.count('-')
        #if Ngap > (0.8 * len(seq)):
        #    continue
        adict[defline] = seq
    return adict


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_distance(name, name2, seqDict):
    global dist_records
    bad_chars = set(['-', 'X', '*'])

    key = [name, name2]
    key.sort()
    key = ':'.join(key)

    try:
        dist = dist_records[key]
    except KeyError:
        seq1 = seqDict[name]
        seq2 = seqDict[name2]
        dist_list = [seq1[i] != seq2[i] for i in range(len(seq1))
                if (seq1[i] not in bad_chars) and (seq2[i] not in bad_chars)]
        if len(dist_list) == 0: #too many gaps or too many similar sequence
            dist = 0
        else:
            dist = float(sum(dist_list)) / len(dist_list)
        dist_records[key] = dist
    return dist


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()


    seqDict = read_fa(fn)
    lengths = [len(seqDict[i]) for i in seqDict]
    if len(set(lengths)) != 1:
        print('alignments are not in the same length! below is the stat:')
        for i in seqDict:
            print(i, len(seqDict[i]))
        sys.exit()


    keys = list(seqDict.keys())

    info = [str(len(seqDict))]
    for name in keys:
        line = [name]
        for name2 in keys:
            if name2 == name:
                line.append(0.0)
            else:
                dist = compute_distance(name, name2, seqDict)
                line.append(dist)
        info.append('\t'.join(map(str, line)))

    info.append('')
    cmn.write_lines(info, fn + '.dist')




