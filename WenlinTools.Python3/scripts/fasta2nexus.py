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

    new = ['#NEXUS\nbegin data;']
    new.append('dimensions ntax=%s nchar=%s;' % (len(seqDict), length))
    new.append('format datatype=DNA interleave=no gap=-;')
    new.append('matrix')
    for name in seqDict:
        new.append('%s        %s' % (name, seqDict[name]))

    new.append(';\nend;\n')

    dn = cmn.lastName(fn) + '.nexus'
    cmn.write_lines(new, dn)


