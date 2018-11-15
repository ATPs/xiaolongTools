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



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py combine.fa.hetero2", file=sys.stderr)
        sys.exit()

    lines = cmn.file2lines(fn)

    header = lines[0].split()[2:]
    #cmn.write_lines(header, 'header_names')

    #new = [' '.join(['sp%s' % i for i in xrange(len(header))])]
    new = [' '.join(header)]
    for line in lines[1:]:
        firstChar = ''
        items = line.split()[2:]
        newline = []
        for item in items:
            chars = item.split()
            if firstChar == '':
                firstChar = chars[0]

            if len(chars) == 1:
                if chars[0] == '-':
                    newline.append('0,0')
                elif firstChar == chars[0]:
                    newline.append('1,0')
                else:
                    newline.append('0,1')

            else:#have both char
                newline.append('1,1')
        new.append(' '.join(newline))

    dn = fn + '.tmix'
    cmn.write_lines(new, dn)

    cmn.run('gzip %s' % dn)



