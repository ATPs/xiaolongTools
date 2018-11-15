#!/usr/bin/env python

#function: tool used to profile a command
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)
import os
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def reformat(info, argvs):
    count = 0
    for argv in argvs:
        count += 1
        label = 'sys.argv[%s]' % count
        info = info.replace(label, '"%s"' % argv)

    info = info.replace("if __name__=='__main__':", "def profiler_main():")

    info += "\n\nimport cProfile\n\n\ncProfile.run('profiler_main()')\n\n"

    return info



if __name__=='__main__':
    #options=parse_options()
    try:
        cmd=sys.argv[1]
    except:
        print("Usage: *.py 'seq2ref.py 254780193'", file=sys.stderr)
        print("the command must contain full python to read it", file=sys.stderr)
        sys.exit()


    import cmn

    argvs = cmd.split()


    info = cmn.txt_read(argvs[0])

    if "__name__=='__main__'" not in info:
        print("program doesn't contain the line: __name__=='__main__'", file=sys.stderr)
        print("exit! do nothing", file=sys.stderr)
        sys.exit()


    #reformat to make it workable for profiler
    info = reformat(info, argvs[1:])

    dn = 'profile_%s' % argvs[0]
    cmn.write_file(info, dn)
    report = cmn.cmd2info('python %s' % dn)

    dn = '%s_report' % argvs[0]
    cmn.write_file(report, dn)
    print('results in %s' % dn)

    dn2 = "%s_sorted" % dn
    cmd = 'cat %s| sort -r -nk4 > %s' % (dn, dn2)
    os.system(cmd)
    print('sorted result by the accumuated time is in %s' % dn2)


