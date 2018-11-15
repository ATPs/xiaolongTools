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
badchars = ['.']

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        #defline, seq = each.strip().split()
        lines = each.strip().split('\n')
        defline = name2ID(lines[0])
        seq = ''.join(lines[1:]).replace('N', '-')
        adict[defline] = seq
    return adict, len(seq)


def name2ID(name):
    ID = name.strip().split('.Ler')[0].split('_')[0]
    return ID

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fnlist, outdir = sys.argv[1:]
    except:
        print("Usage: *.py falist outdirName", file=sys.stderr)
        sys.exit()

    outlabel = cmn.lastName(fnlist)
    #scan through to see the total IDlist
    IDs = set([])
    fns = cmn.file2lines(fnlist)
    for fn in fns:
        with open(fn) as fp:
            for line in fp:
                if line[0] == '>':
                    ID = name2ID(line[1:].strip())
                    IDs.add(ID)

    #read in sequence and partition
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

    outdir = outdir + '_pfind'
    cmn.mkdir(outdir)
    os.chdir(outdir)

    #output the phylip format file
    seqDict = {ID: ''.join(final[ID]) for ID in final}
    length = len(seqDict[ID])

    new = ['%s\t%s' % (len(seqDict), length)]
    for name in seqDict:
        new.append('%s        %s' % (name, seqDict[name]))

    dn = outlabel + '.phylip'
    cmn.write_lines(new, dn)

    #write out the partition file
    ftemplate = '/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/partition_finder.cfg.template'
    fcfg = 'partition_finder.cfg'
    info = cmn.txt_read(ftemplate)
    info = info.replace('[input_phylip]', dn)
    ##Gene3_pos3 = 1452-2208\3;
    print('assuming all are protein coding genes')
    blocks = []
    for name,i,j in setList:
        for pos in [0, 1, 2]:
            pLabel = '%s_%s' % (name, pos+1)
            for char in badchars:
                pLabel = pLabel.replace(char, '_')
            line = '%s = %s-%s\\3;\n' % (pLabel, i+pos, j)
            blocks.append(line)
    info = info.replace('[data_block_input]', ''.join(blocks))
    cmn.write_file(info, fcfg)

    print('please open the config file to make sure it is correct:')
    print('%s/%s' % (outdir, fcfg))

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

    dn = outlabel + '.nexus'
    cmn.write_lines(new, dn)


