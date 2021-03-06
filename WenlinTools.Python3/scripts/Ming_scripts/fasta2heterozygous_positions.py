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
def chars2line(pos, chars):
    pos = str(pos)
    line = [pos, pos]
    for i in range(0, len(chars), 2):
        a, b = chars[i:i+2]
        if a == b:
            line.append(a)
        else:
            line.append('%s,%s' % (a, b))
    return '\t'.join(line)


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

    adict = {}
    with open(fn) as fp:
        for line in fp:
            if line[0] == '>':
                name = line[1:].split('_')[0]
            else:
                seq = line.strip()
                try:
                    adict[name].append(seq)
                except KeyError:
                    adict[name] = [seq]

    names = list(adict.keys())
    header = ['fake', 'position'] + names
    header = '\t'.join(header)

    seqs = []
    hetero = [[header], [header], [header], [header]]
    for name in names:
        seqs += adict[name]

    for i in range(len(seqs[0])):
        chars = [seq[i] for seq in seqs]
        check_chars = [char for char in chars if (char != '-') and (char != 'N')]
        N = len(set(check_chars))
        if N == 0:
            #skip the all gapped positions
            continue
        else:
            line = chars2line(i, chars)
            hetero[N-1].append(line)


    for i, each in enumerate(hetero):
        dn = cmn.lastName(fn) + '.hetero%s' % (i + 1)
        cmn.write_lines(each, dn)



