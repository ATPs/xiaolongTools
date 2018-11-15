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

rdict = {
        'A': 'T',
        'T': 'A',
        'G': 'C',
        'C': 'G',
        'N': 'N',
        '-': 'N',
        '*': 'N'
        }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        read, fn, direction = sys.argv[1:4]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    fished_reads = []
    cmn.mkdir('grep_out')
    if True:
        with open(fn) as fp:
            for lineN, line in enumerate(fp):
                if lineN % 4 != 1:#only take the sequence
                    continue
                line = line.strip()
                if direction == 'backward':
                    line = line[::-1]

                #find match forward and + strand
                i1 = line.find(read)
                if i1 != -1:
                    fished_reads.append(line[i1:])


                #in the reverse strand
                try:
                    rline = ''.join([rdict[i] for i in line[::-1]])
                except:
                    print('cannot reverse the read: %s' % line)
                    continue

                line = rline

                i3 = line.find(read)
                if i3!= -1:
                    fished_reads.append(line[i3:])

    cmn.write_lines(fished_reads, 'grep_out/%s.grep' % cmn.lastName(fn))


