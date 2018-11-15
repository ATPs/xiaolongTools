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
        print("Usage: *.py *.fa", file=sys.stderr)
        sys.exit()



    cmds = ['export PATH=$PATH:/home2/wli/local/EIG6.1.1/bin;']
    cmds.append('module add gsl/1.15;')
    cmds.append('module add openblas/intel/0.2.14;')

    #transform_cmd = 'fasta2EIGinput.py %s;' % fn
    #cmds.append(transform_cmd)

    #fn = fn + '_loose'

    maincmd = '/home2/wli/local/EIG6.1.1/bin/smartpca.pl -i %s.geno -a %s.snp -b %s.ind' % (fn, fn, fn)

    dnlabel = cmn.lastName(fn).replace('.fa', '')
    maincmd += ' -o %s.pca -p %s.plot -e %s.eval -l %s.log' % (dnlabel, dnlabel, dnlabel, dnlabel)

    cmds.append(maincmd)

    cmn.run('\n'.join(cmds))




