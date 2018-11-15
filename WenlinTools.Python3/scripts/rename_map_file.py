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
        reflabel = sys.argv[2]
    except:
        print("Usage: *.py maplist reflabel", file=sys.stderr)
        sys.exit()


    for fn in cmn.file2lines(fn):
        items = fn.split('/')
        if len(items) == 1:
            wdir = '.'
        else:
            wdir = '/'.join(items[:-1])
        lastName = items[-1]
        sample = lastName.split('_')[0]
        newName = '%s_%s_snp_step2.map' % (sample, reflabel)
        cmd = 'cd %s; mv %s %s' % (wdir, lastName, newName)
        #print 'change %s to %s' % (lastName, newName)
        #cmn.run(cmd)
        print(cmd)





