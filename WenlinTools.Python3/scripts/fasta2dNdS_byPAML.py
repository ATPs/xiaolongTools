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
import ete3
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def prune_tree(ftree, fseq):
    t = ete3.Tree(ftree)
    IDlist = cmn.cmd2lines('grep ">" %s|cut -d ">" -f 2' % fseq)
    t.prune(IDlist)
    dn = 'prune_tree.tre'
    cmn.write_file(t.write(format=1), dn)
    return dn


if __name__=='__main__':
    #options=parse_options()
    try:
        fn, ftree, outlabel = sys.argv[1:]
        ftree = os.path.abspath(ftree)
        fn = os.path.abspath(fn)
    except:
        print("Usage: *.py fasta ftree outlabel", file=sys.stderr)
        sys.exit()

    #prepare a temp directory to run the code
    tmpdir = 'tmp_%s' % outlabel
    cmn.mkdir(tmpdir)
    os.chdir(tmpdir)

    #make the nuc input
    fnuc = cmn.lastName(fn) + '.nuc'
    cmd = '/project/biophysics/Nick_lab/wli/sequencing/scripts/fasta2pamlCodonAlignment.py %s' % fn
    cmn.run(cmd)

    #prune the tree by the shown taxons
    ftree = prune_tree(ftree, fn)

    #make the configuration
    ftemplate = '/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/CODEML-M0.ctl.template'
    info = cmn.txt_read(ftemplate)
    info = info.replace('[input]', fnuc)
    info = info.replace('[ftree]', ftree)
    info = info.replace('[outlabel]', outlabel)

    fctl = 'codeM.ctl'
    cmn.write_file(info, fctl)

    #run the code
    cmd = '/home2/wli/local/paml4.9g/src/codeml codeM.ctl > codeml.log'
    cmn.run(cmd)

    #clean up
    fout = outlabel + '.output'
    cmd = 'cp %s ../' % fout
    cmn.run(cmd)

    os.chdir('..')
    cmn.run('rm -r %s' % tmpdir)




