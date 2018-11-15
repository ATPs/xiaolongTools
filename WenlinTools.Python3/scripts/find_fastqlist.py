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

def nonGap_char(fn):
    seq = ''.join([line.strip() for line in cmn.file2lines(fn)
        if not line[0] == '>'])
    seq = seq.replace('N', '-')
    N = len(seq) - seq.count('-')
    return N


def tell_best_fa(alist):
    print('warnning: telling best fasta from the following sequences')
    print('\n'.join(alist))
    bestOne = max(alist, key=lambda x: nonGap_char(x))
    print('the best one is %s' % bestOne)
    return bestOne
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def count_ass_appearance(alist):
    #LEP14801_3318_gapClosed_withMito_cut1000_snp_step2_m2s.fa
    adict = {}
    for fn in alist:
        items = cmn.lastName(fn).split('_')
        ass = '_'.join(items[1:-3]).replace('_withMito', '')
        try:
            adict[ass] += 1
        except KeyError:
            adict[ass] = 1
    return adict


def tell_best_falist(alist):
    #find pairs
    R1, R2, single = [], [], []
    nonAlea_list = [each for each in alist
            if '/archive/butterfly' not in each]
    if len(nonAlea_list) != 0:
        return nonAlea_list
    return alist


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py NickList requred_keywords", file=sys.stderr)
        sys.exit()

    words = sys.argv[2:]

    IDs = set([line.strip().split()[0] for line in cmn.file2lines(fn)
        if line.strip() != ''])
    #TODO: find common references

    wdirs = [
            '/project/biophysics/Nick_lab/mtang/archive/fastq_libs',
            '/project/biophysics/Nick_lab/mtang/archive/reference_fastq/*',
            '/project/biophysics/Nick_lab/mtang/preprocessing/*',
            '/archive/biophysics/Nick_lab/mtang/preprocessing/*',
            '/project/biophysics/Nick_lab/jzhang/wenlin/mapping/*/*',
            '/work/biophysics/mtang/preprocessing/*',
            ]

    missing = []
    falist = []
    required = ''
    if len(words) != 0:
        required = '|' + '|'.join(['grep %s' % word for word in words])

    alea_list = cmn.cmd2lines('ssh butterfly@toxea.swmed.edu "ls /archive/butterfly/ready_fastq/*q"')
    alea_list += cmn.cmd2lines('ssh butterfly@toxea.swmed.edu "ls /archive/butterfly/reference_fastq/*/*q"')
    #alea_list += cmn.cmd2lines('ssh wenlin@alea.swmed.edu "ls /home/jshen/h7/*/*/*q"')
    #alea_list += cmn.cmd2lines('ssh wenlin@alea.swmed.edu "ls /home/jshen/h8/*/*/*q"')


    for ID in IDs:
        taken = []
        for wdir in wdirs:
            cmd = 'ls %s/%s_*q  2> /dev/null| grep -v MITO %s' % (wdir, ID, required)
            taken += cmn.cmd2lines(cmd)

        taken += [each for each in alea_list
                if cmn.lastName(each).split('_')[0] == ID]

        if len(taken) < 2:
            missing.append(ID)
            continue
        else:
            taken = tell_best_falist(taken)

        falist += taken

    #final filtering to remove cases where finding 1510 would get 15101C03
    tmp = []
    checklist = list(IDs)
    for line in falist:
        sp = cmn.lastName(line).split('_')[0]
        if sp in IDs:
            tmp.append(line)
            try:
                checklist.remove(sp)
            except:
                pass

    #print falist
    falist = tmp
    #print falist
    missing = set(checklist) | set(missing)

    alea_files = [each for each in falist
            if '/archive/butterfly/' in each or ('jshen/h' in each)]

    if len(alea_files) != 0:
        print('the following files need to transfer from /archive server')
        print('\n'.join(alea_files))

    cmn.write_lines(falist, 'attempt.fastqlist')

    if len(missing) != 0:
        print('ATTENTION! the following ID missing sequence!')
        print('\n'.join(missing))
        cmn.write_lines(missing, 'missingIDs')

