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
        print("Usage: *.py falist_file", file=sys.stderr)
        sys.exit()

    fns = cmn.getid(fn)
    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/RAxML_tree.job')
    count=0
    for fn in fns:
        count += 1
        label = cmn.lastName(fn).replace('.fa', '')
        info = template.replace('[FN]', fn)
        info = info.replace('[outlabel]', label)
        dn = 'BTtree%s.job' % count
        cmn.write_file(info, dn)




