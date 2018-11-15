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
        sp=sys.argv[1]
    except:
        print("Usage: *.py spID", file=sys.stderr)
        sys.exit()


    fns = ['rescued_ratio_mis1.txt', 'good_read_assembled.txt', 'final_comparison.table']
    fmis = cmn.cmd2lines('ls rescued_read_assembled_mis1*.txt')[0]
    fns.append(fmis)

    for fn in fns:
        cmd = 'chmod a+w %s' % fn
        cmn.run(cmd)

    cmd = "ssh wenlin@morpho.swmed.edu 'rm /data/www/wenlin/html/transfer/barcode_lineup_files/%s_rescued_read_assembled_mis1*.txt'" % sp
    cmn.run(cmd)

    cmd = 'rm /project/biophysics/Nick_lab/wli/archive/BWA_barcodes/lineup_files/%s_rescued_read_assembled_mis1*.txt' % sp
    cmn.run(cmd)

    ddirs = ['/project/biophysics/Nick_lab/wli/archive/BWA_barcodes/lineup_files',
            'wenlin@morpho.swmed.edu:/data/www/wenlin/html/transfer/barcode_lineup_files/other_data']

    cmd = 'rsync -av %s wenlin@morpho.swmed.edu:/data/www/wenlin/html/transfer/barcode_lineup_files/%s_%s' % (fmis, sp, fmis)
    cmn.run(cmd)

    for ddir in ddirs:
        for fn in fns:
            if 'rescued_read_assembled_mis1' in fn and '/other_data' in ddir:
                continue

            cmd = 'rsync -av %s %s/%s_%s' % (fn, ddir, sp, cmn.lastName(fn))
            cmn.run(cmd)


