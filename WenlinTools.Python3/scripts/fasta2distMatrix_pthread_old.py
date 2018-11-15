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


    keys = list(seqDict.keys())

    result_list = []
    for name in keys:
        process = pool.apply_async(compute_one_line, [name, keys, seqDict])
        result_list.append(process)


    rdict = {}
    for process in result_list:
        [name, line] = process.get()
        rdict[name] = line


    info = [str(len(seqDict))]
    for name in keys:
        info.append(rdict[name])

    info.append('')
    cmn.write_lines(info, fn + '.dist')


