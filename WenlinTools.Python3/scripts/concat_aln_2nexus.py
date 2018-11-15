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
        defline = lines[0].strip().split('_')[0]
        seq = ''.join(lines[1:]).replace('N', '-')
        adict[defline] = seq
    return adict, len(seq)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fnlist=sys.argv[1]
    except:
        print("Usage: *.py falist", file=sys.stderr)
        sys.exit()

    IDs = set([])
    fns = cmn.file2lines(fnlist)
    for fn in fns:
        with open(fn) as fp:
            for line in fp:
                if line[0] == '>':
                    ID = line[1:].strip().split('_')[0]
                    IDs.add(ID)

    shift = 0
    final = {}
    setList = []
    for fn in fns:
        totalIDs = set(IDs)
        namelabel = '.'.join(cmn.lastName(fn).split('.')[:-1])
        seqDict, seqLength = read_fa(fn)
        setList.append((namelabel, shift + 1, shift + seqLength))
        shift += seqLength
        for ID in seqDict:
            line = seqDict[ID]
            if True:
                try:
                    final[ID].append(line)
                except KeyError:
                    final[ID] = [line]

        missingIDs = totalIDs - set(seqDict.keys())
        #print totalIDs, seqDict.keys()
        #fill in missing IDs
        for ID in missingIDs:
            line = '-' * seqLength
            try:
                final[ID].append(line)
            except KeyError:
                final[ID] = [line]

        seqDict, length = read_fa(fn)

    new = ['#NEXUS\nbegin data;']
    new.append('dimensions ntax=%s nchar=%s;' % (len(final), shift))
    new.append('format datatype=DNA interleave=no gap=-;')
    new.append('matrix')
    for name in final:
        new.append('%s        %s' % (name, ''.join(final[name])))

    new.append(';\nend;\n')

    new.append('begin sets;\n')
    geneInfo = []
    count = 0
    for name,i,j in setList:
        new.append('charset %s = %s-%s' % (name, i, j))
        geneInfo.append('s%s:%s' % (count, name))
        count += 1

    new.append('charpartition byGene = %s;' % (', '.join(geneInfo)))
    new.append('\nend;\n')

    dn = 'concat.nxs'
    cmn.write_lines(new, dn)


