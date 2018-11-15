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

    final = {}
    for fn in fns:
        for line in cmn.file2lines(fn):
            if line.strip() == '':
                continue

            if line[0] == '>':
                ID = line.strip()
            else:
                try:
                    final[ID].append(line)
                except KeyError:
                    final[ID] = [line]
    dn = 'concat.fasta'
    fasta = ['%s\n%s\n' % (name, ''.join(final[name]))
            for name in final]
    cmn.write_file(''.join(fasta), dn)



