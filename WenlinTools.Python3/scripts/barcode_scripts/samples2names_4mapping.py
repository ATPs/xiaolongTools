#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/work/00412/mtang/sequencing/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
#import pandas as pd
#import datetime
from fullname_lib import get_names_4barcode

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        #fn, f_table = sys.argv[1:3]
        fn = sys.argv[1]
    except:
        print("Usage: *.py RAxML_bestTree.noGap", file=sys.stderr)
        sys.exit()

    nameDict = get_names_4barcode()

    info = []
    missing = []
    lines = cmn.getid(fn)

    for line in lines:
        sp = line.strip().split()[0]
        #line = '%s\t%s\n' % (line, nameDict[sp])
        try:
            info.append(nameDict[sp].replace('"', ''))
        except KeyError:
            missing.append(sp)

    info.append('')

    info = '\n'.join(info)

    cmn.write_file(info, 'sampleInfo')

    print('caution! the following IDs are missing:')
    print('\n'.join(missing))

