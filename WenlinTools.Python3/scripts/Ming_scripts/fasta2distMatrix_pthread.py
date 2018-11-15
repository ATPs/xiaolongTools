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
from multiprocessing import Pool
dist_records = {}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        #defline, seq = each.strip().split()
        #seq = seq.strip()
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
    bad_chars = set(['-', 'N'])

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
        dist = float(sum(dist_list)) / len(dist_list)
        dist_records[key] = dist
    return dist


def compute_one_line(name, keys, seqDict):
    if True:
        line = [name]
        for name2 in keys:
            if name2 == name:
                line.append(0.0)
            else:
                dist = compute_distance(name, name2, seqDict)
                line.append(dist)

    return [name, '\t'.join(map(str, line))]


def compute_distance_for_piece(plabel, seq1, seq2):
    bad_chars = set(['-', 'N'])
    length = len(seq1)

    check = [seq1[i] != seq2[i] for i in range(length)
            if (seq1[i] not in bad_chars) and (seq2[i] not in bad_chars)]

    return plabel, sum(check), len(check)


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, cpu = sys.argv[1:3]
        cpu = int(cpu)
    except:
        print("Usage: *.py *.fa cpu", file=sys.stderr)
        sys.exit()

    pool = Pool(processes=cpu)              # start 4 worker processes

    seqDict = read_fa(fn)
    lengths = [len(seqDict[i]) for i in seqDict]
    if len(set(lengths)) != 1:
        print('alignments are not in the same length! below is the stat:')
        for i in seqDict:
            print(i, len(seqDict[i]))
        sys.exit()

    length = lengths[0]

    keys = list(seqDict.keys())

    #make the distance dict first using pthreads
    result_list = []
    rdict = {}
    bunch = 10000 # compute them for ever 100 letters
    for i, name1 in enumerate(keys):
        seq1 = seqDict[name1]
        for name2 in keys[i+1:]:
            seq2 = seqDict[name2]
            for ii in range(0, length, bunch):
                piece1 = seq1[ii: ii + bunch]
                piece2 = seq2[ii: ii + bunch]
                plabel = (name1, name2)
                process = pool.apply_async(compute_distance_for_piece, [plabel, piece1, piece2])
                result_list.append(process)
    
    pool.close()
    pool.join()

    for process in result_list:
        [plabel, Ndiff, Ntotal] = process.get()
        if plabel not in rdict:
            rdict[plabel] = [0, 0]

        rdict[plabel][0] += Ndiff
        rdict[plabel][1] += Ntotal

    dist_dict = {key: float(rdict[key][0])/ rdict[key][1]
            for key in rdict}


    info = [str(len(seqDict))]
    for i, name in enumerate(keys):
        line = [name]
        for j, name2 in enumerate(keys):
            if name2 == name:
                line.append(0.0)
            else:
                if j > i:
                    plabel = (name, name2)
                elif j < i:
                    plabel = (name2, name)
                dist = dist_dict[plabel]
                line.append(dist)
        info.append('\t'.join(map(str, line)))

    info.append('')
    cmn.write_lines(info, fn + '.dist')




