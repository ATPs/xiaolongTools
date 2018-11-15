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
        fn=sys.argv[1]
        fn = os.path.abspath(fn)
    except:
        print("Usage: *.py *.fa [ignore]", file=sys.stderr)
        sys.exit()

    ignore_check = False
    if len(sys.argv) > 2:
        if sys.argv[2] == 'ignore':
            ignore_check = True

    outlabel = cmn.lastName(fn).replace('.fa', '')
    wdir = '%s_tmp' % cmn.lastName(fn)
    cmn.mkdir(wdir)
    os.chdir(wdir)

    fphy = cmn.lastName(fn) + '.phylip'
    fname = cmn.lastName(fn) + '.phylipNames.dict.pkl'

    fchecks = ['outfile', 'dist.Tree']
    isbad = False
    for fcheck in fchecks:
        if not ignore_check and cmn.filexist(fcheck):
            print('Erorr: file %s exists! running pipeline would overwrite the files, please either delete it or move it to another place' % fcheck)
            isbad = True
    if isbad:
        sys.exit()

    cmd = 'source /home2/wli/.bash_profile;/project/biophysics/Nick_lab/wli/sequencing/scripts/fasta2phylip4dnadist.py %s' % fn
    cmn.run(cmd)

    dnadistInfo = '%s\nY\n' % fphy
    cmn.write_file(dnadistInfo, 'input.dnadist')

    cmd = 'rm outfile 2> /dev/null;/home2/wli/local/phylip-3.696/exe/dnadist < input.dnadist > dnadist.log'
    #print cmd
    cmn.run(cmd)

    bionjInfo = 'outfile\ndist.Tree\n'
    cmn.write_file(bionjInfo, 'input.bionj')
    cmd = 'source /home2/wli/.bash_profile; BIONJ_linux < input.bionj > bionj.log'
    cmn.run(cmd)

    cmd = 'source /home2/wli/.bash_profile; /home2/wli/anaconda/bin/python /project/biophysics/Nick_lab/wli/sequencing/scripts/rename_phylipnaming.py dist.Tree %s' % fname
    cmn.run(cmd)


    dn = '%s_dnadistBIONJ.tre' % outlabel
    dmatrix = '%s_dnadist.matrix' % outlabel
    cmd = 'mv outfile ../%s; mv dist.Tree.mapnamed ../%s; mv %s ..' % (dmatrix, dn, fname)
    cmn.run(cmd)

    os.chdir('..')
    #cmn.run('rm -r ' +  wdir)
