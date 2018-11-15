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
        IDlist = cmn.file2lines(sys.argv[2])
    except:
        print("Usage: *.py ind_file takenIDs", file=sys.stderr)
        sys.exit()


    new = []
    for line in cmn.file2lines(fn):
        items = line.strip().split()
        Id = items[0]
        if Id in IDlist:
            items[-1] = 'case1'
        else:
            items[-1] = 'Ignore'

        new.append('\t'.join(items))

    new.append('')
    cmn.write_lines(new, fn)




