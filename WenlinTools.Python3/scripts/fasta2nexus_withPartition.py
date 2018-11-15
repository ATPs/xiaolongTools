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

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        #defline, seq = each.strip().split()
        lines = each.strip().split('\n')
        defline = lines[0]
        seq = ''.join(lines[1:])
        seq = seq.upper()
        adict[defline] = seq
    return adict, len(seq)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def check_alignment_length(fn):
    taken = []
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                if len(taken) == 0:
                    continue
                else:
                    break
            else:
                taken.append(line.strip())
    length = len(''.join(taken))
    return length

def fn2name(fn):
    try:
        name = fn.split('/')[-1]
    except:
        name = cmn.lastName(fn)
    return name


def parse_falist(fn):
    lines = cmn.file2lines(fn)
    rlist = []
    lastN = 0
    for fn in lines:
        if fn[0] == '#':
            continue

        length = check_alignment_length(fn)
        name = fn2name(fn)
        i = lastN + 1
        j = lastN + length
        lastN = j
        rlist.append((name, i, j))

    return rlist



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, f_falist=sys.argv[1:]
    except:
        print("Usage: *.py fa falist.ToConcated", file=sys.stderr)
        print("the falist.ToConcated is used to check the partition ordering", file=sys.stderr)
        sys.exit()

    seqDict, length = read_fa(fn)
    setList = parse_falist(f_falist)

    new = ['#NEXUS\nbegin data;']
    new.append('dimensions ntax=%s nchar=%s;' % (len(seqDict), length))
    new.append('format datatype=DNA interleave=no gap=-;')
    new.append('matrix')
    for name in seqDict:
        new.append('%s        %s' % (name, seqDict[name]))

    new.append(';\nend;\n')

    new.append('begin sets;\n')
    geneInfo = []
    count = 0
    for name,i,j in setList:
        new.append('charset %s = %s-%s;' % (name, i, j))
        geneInfo.append('s%s:%s' % (count, name))
        count += 1

    new.append('charpartition byGene = %s;' % (', '.join(geneInfo)))
    new.append('\nend;\n')

    dn = cmn.lastName(fn) + '.nexus'
    cmn.write_lines(new, dn)


