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
mash_file_dict = {}

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



def get_mash_file(name, seq):
    global mash_file_dict, cpu
    try:
        fn = mash_file_dict[name]
    except KeyError:
        fn = '/tmp/%s' % name
        seq = ''.join(seq).replace('-', '').replace('N', '')
        fasta = '>%s\n%s\n' % (name, seq)
        cmn.write_file(fasta, fn)
        cmd = '/home2/wli/local/mash-Linux64-v1.1.1/mash sketch -n -p %s %s' % (cpu, fn)
        cmn.run(cmd)
        dn = fn + '.msh'
        mash_file_dict[name] = dn
        fn = dn
    return fn


def compute_mash_distance(f1, f2):
    global cpu
    dn = '/tmp/%s-%s' % (cmn.lastName(f1), cmn.lastName(f2))
    cmd = '/home2/wli/local/mash-Linux64-v1.1.1/mash dist -p %s %s %s > %s' % (cpu, f1, f2, dn)
    cmn.run(cmd)
    #print 'cmd:', cmd
    dist = cmn.txt_read(dn).strip().split()[2]
    return dist

def mash_distance(name1, name2, seq1, seq2):
    f1 = get_mash_file(name1, seq1)
    f2 = get_mash_file(name2, seq2)
    dist = compute_mash_distance(f1, f2)
    return dist
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def compute_distance_byMash(name, name2, seqDict):
    global dist_records
    #bad_chars = set(['-', 'N'])

    key = [name, name2]
    key.sort()
    key = ':'.join(key)

    try:
        dist = dist_records[key]
    except KeyError:
        seq1 = seqDict[name]
        seq2 = seqDict[name2]
        dist = mash_distance(name, name2, seq1, seq2)
        #dist_list = [seq1[i] != seq2[i] for i in xrange(len(seq1))
        #        if (seq1[i] not in bad_chars) and (seq2[i] not in bad_chars)]
        #if len(dist_list) == 0: #too many gaps or too many similar sequence
        #    dist = 0
        #else:
        #    dist = float(sum(dist_list)) / len(dist_list)
        dist_records[key] = dist
    return dist


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
        cpu = int(sys.argv[2])
    except:
        print("Usage: *.py *.fa cpu", file=sys.stderr)
        sys.exit()


    seqDict = read_fa(fn)
    lengths = [len(seqDict[i]) for i in seqDict]
    if len(set(lengths)) != 1:
        print('Warnning: alignments are not in the same length! below is the stat:')
        for i in seqDict:
            print(i, len(seqDict[i]))


    keys = list(seqDict.keys())

    info = [str(len(seqDict))]
    for name in keys:
        line = [name]
        for name2 in keys:
            if name2 == name:
                line.append(0.0)
            else:
                dist = compute_distance_byMash(name, name2, seqDict)
                line.append(dist)
        info.append('\t'.join(map(str, line)))

    info.append('')
    cmn.write_lines(info, fn + '.mash.dist')




