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

    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/beast_template.xml')

    info = []
    for name in seqDict:
        info.append('<sequence taxon="%s">%s</sequence>' % (name, seqDict[name]))

    info.append('')
    new = template.replace('[WLdata]', '\n'.join(info))
    new = new.replace('[WLlabel]', cmn.lastName(fn))

    dn = cmn.lastName(fn) + '.xml'
    cmn.write_file(new, dn)


