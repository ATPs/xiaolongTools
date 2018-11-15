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
        print("Usage: *.py fa", file=sys.stderr)
        print('this code would use 12 CPU in parallel', file=sys.stderr)
        sys.exit()

    dn = cmn.lastName(fn) + '.modelCheck'
    cmd = 'java -jar /home2/wli/local/jmodeltest-2.1.10/jModelTest.jar '
    cmd += '-d %s -AIC -f -g 6 > %s' % (fn, dn)
    cmn.run(cmd)




