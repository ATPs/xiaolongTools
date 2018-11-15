#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

import cmn

if __name__=='__main__':
    #options=parse_options()
    try:
        script=sys.argv[1]
        fn=sys.argv[2]
        N=int(sys.argv[3])
    except:
        print("Usage: *.py main_script parameter_file parts [tmpdir:tmpx]", file=sys.stderr)
        sys.exit()

    try:
        tmpdir = sys.argv[4] + 'tcmds'
    except:
        tmpdir = 'tmpx'
    paras = cmn.file2lines(fn)
    length = len(paras)
    bundle = length / N

    cmn.mkdir(tmpdir)

    cmds = []
    for i in range(N):
        lines = paras[i*bundle:(i+1)*bundle]
        if i == N - 1:
            lines += paras[(i+1)*bundle:]
        dn = '%s/tmp_%s' % (tmpdir, i)
        cmn.write_file('\n'.join(lines), dn)
        cmds.append('%s %s' % (script, dn))


    cmn.write_file('\n'.join(cmds), '%s_parsed_cmds' % tmpdir)

