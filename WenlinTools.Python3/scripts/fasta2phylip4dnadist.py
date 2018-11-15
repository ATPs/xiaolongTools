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
        seq = ''.join(lines[1:]).replace('N', '-')
        adict[defline] = seq
    return adict, len(seq)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    seqDict, length = read_fa(fn)
    nameDict = {}
    count = 0

    new = ['%s\t%s' % (len(seqDict), length)]
    for name in seqDict:
        count += 1
        newName = 'ID%s' % (count)
        nameDict[newName] = name
        newName = '{:<10}'.format(newName)
        new.append('%s%s' % (newName, seqDict[name]))


    dn = cmn.lastName(fn) + '.phylip'
    cmn.write_lines(new, dn)

    dn = cmn.lastName(fn) + '.phylipNames.dict.pkl'
    cmn.pickle_write(nameDict, dn)

