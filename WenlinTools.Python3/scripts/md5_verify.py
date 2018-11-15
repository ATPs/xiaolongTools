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
import hashlib

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn1 = sys.argv[1]
        fn2 = sys.argv[2]
    except:
        print("Usage: *.py fn fn.md5", file=sys.stderr)
        sys.exit()


    md5_string = cmn.txt_read(fn2).split()[0]

    info = cmn.txt_read(fn1)

    m = hashlib.md5()

    m.update(info)

    fn_string = m.hexdigest()

    if fn_string != md5_string:
        print('Error! md5 doesn\'t match for %s' % fn1)
    else:
        print('success! md5 matches for %s' % fn1)



