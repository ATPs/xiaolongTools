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



def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        ID = name2ID(defline)
        seq = ''.join(lines[1:])
        adict[ID] = seq
    return adict, len(seq)

def name2ID(name):
    return name.strip()
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py filelist", file=sys.stderr)
        sys.exit()


    fns = cmn.getid(fn)
    #check once to get all IDs
    IDs = set([])
    for fn in fns:
        with open(fn) as fp:
            for line in fp:
                if line[0] == '>':
                    ID = name2ID(line[1:])
                    IDs.add(ID)

    final = {}
    for fn in fns:
        print(fn)
        totalIDs = set(IDs)
        seqDict, seqLength = read_fa(fn)
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

    dn = 'concat.fasta'
    fasta = ['>%s\n%s\n' % (name, ''.join(final[name]))
            for name in final]
    cmn.write_file(''.join(fasta), dn)



