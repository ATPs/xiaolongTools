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
import os


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn, K = sys.argv[1:3]
        K = int(K)
    except:
        print("Usage: *.py mainparams K", file=sys.stderr)
        print("mainparams were generated using GUI", file=sys.stderr)
        sys.exit()


    config = cmn.file2lines(fn)
    rep = 4
    newConfig = []
    cwd = os.getcwd()

    for line in config:
        if '#define OUTFILE' in line:
            newConfig.append('#define OUTFILE structure.output')
        elif '#define MAXPOPS' in line:
            newConfig.append('#define MAXPOPS %s' % K)
        else:
            newConfig.append(line)

    newConfig.append('')

    cmds = []
    for eachtime in range(rep):
        outdir = 'structureK%s/r%s' % (K, eachtime)
        cmn.mkdir(outdir)
        os.chdir(outdir)

        dn = 'mainparams'
        cmn.write_lines(newConfig, dn)

        cmd = 'cd %s; touch extraparams; /home2/wli/local/Structure/bin/structure > runtime.log & cd %s' % (outdir, cwd)
        #cmn.run(cmd)

        cmds.append(cmd)
        cmd = 'sleep 55;'
        cmds.append(cmd)

        os.chdir(cwd)

    cmds.append('\nwait\n')
    dn = 'masterK%s.bash' % K
    cmn.write_lines(cmds, dn)
    cmd = 'bash %s' % dn
    cmn.run(cmd)

