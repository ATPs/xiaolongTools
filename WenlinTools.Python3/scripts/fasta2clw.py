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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def read_fa(fa):
    adict = {}
    fastas = cmn.txt_read(fa).split('>')[1:]
    alist = []
    for each in fastas:
        lines = each.strip().split('\n')
        defline = lines[0]
        #defline = defline.replace('-', '_').replace('?', '')
        defline = defline.replace('?', '')
        seq = ''.join(lines[1:])
        alist.append(defline)
        adict[defline] = seq
    return adict, alist




if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()

    #snice this code was used to get barcode for Nick, add a function of checking alignment length

    dn = cmn.lastName(fn).replace('.fasta', '').replace('.fa', '') + '.clw'

    seqDict, orderlist = read_fa(fn)

    maxLength = max(list(map(len, list(seqDict.keys()))))

    strFormat = '{:<%s}' % maxLength

    new = []
    lenDict = {}
    for name in orderlist:
        seq = seqDict[name]
        try:
            lenDict[len(seq)].append(name)
        except:
            lenDict[len(seq)] = [name]
        name = name.replace(' ', '_').replace('\t', '_')
        line = '%s    %s\n' % (strFormat.format(name), seq)
        new.append(line)

    cmn.write_file(''.join(new), dn)


    if len(lenDict) != 1:
        print('Warning: more than one sequence length in the alignment')
        for length in lenDict:
            names = lenDict[length]
            print('\nlength %s: %s\n' % (length, ', '.join(names)))
