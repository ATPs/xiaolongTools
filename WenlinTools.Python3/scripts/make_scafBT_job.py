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

def transfer_alea_files(fnlist):
    global transferDir
    newlist = []
    for fn in fnlist:
        print('transfering %s from archive server ...' % fn)
        cmd = 'rsync -r butterfly@alea.swmed.edu:%s %s' % (fn, transferDir)
        cmn.run(cmd)
        newlist.append('%s/%s' % (transferDir, cmn.lastName(fn)))
    return newlist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        fn = sys.argv[1]
    except:
        print("Usage: *.py", file=sys.stderr)
        sys.exit()

    #check if all the files has contains
    falist = cmn.file2lines(fn)
    bad_falist = [fa for fa in falist
            if not cmn.filexist(fa) and '/archive/butterfly/' not in fa]

    if len(bad_falist) != 0:
        print('Error!')
        print('the following files are errorous:')
        print('\n'.join(bad_falist))
        sys.exit()

    transferDir = 'archiveTransfer'
    cmn.mkdir(transferDir)


    alea_list = [fa for fa in falist
            if '/archive/butterfly' in fa]

    biohpc_list = set(falist) - set(alea_list)

    newlist = transfer_alea_files(alea_list)

    newlist += list(biohpc_list)

    dn = 'new.falist'
    cmn.write_lines(newlist, dn)
    #backup this newlist
    cmn.mkdir('../falist_info')
    dirlabel = os.getcwd().rstrip('/').split('/')[-1]
    backFn = '../falist_info/%s.falist' % dirlabel
    cmd = 'cp %s %s' % (dn, backFn)
    cmn.run(cmd)

    fass = 'all_genomes.fa'

    cmd = 'cat %s| xargs cat > %s' % (dn, fass)
    print('combining fasta into one alignment...')
    cmn.run(cmd)

    #/project/biophysics/Nick_lab/mtang/unbias_SNPs/OK_2016-09-23/step4_postprocessing/map2fasta/6409_Calycopis_cecrops_assembly_V1.1_snp_step2_m2s.fa
    fa = falist[0]
    reflabel = '_'.join(cmn.lastName(fa).split('_')[1:-3])

    fheader = '/work/biophysics/mtang/SNP_calling/indexed_references/%s_scaf.header' % reflabel

    incmd = 'python /project/biophysics/Nick_lab/wli/sequencing/scripts/partition_genomes_scafOrder_4bootstrap_new.py %s %s 100 0.3' % (fass, fheader)

    cmd = 'make_job.py "%s" -p 256GB > scafBT.job' % incmd
    cmn.run(cmd)

    print('please submit job scafBT.job to queue')

    cmd = 'rm -r %s' % transferDir
    cmn.run(cmd)



