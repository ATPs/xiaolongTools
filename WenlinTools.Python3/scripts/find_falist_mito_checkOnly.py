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

def transfer_alea_files(fnlist):
    transferDir = 'alea_transfer'
    cmn.mkdir(transferDir)
    newlist = []
    for fn in fnlist:
        print('transfering %s from archive server ...' % fn)
        cmd = 'rsync -r butterfly@alea.swmed.edu:%s %s' % (fn, transferDir)
        cmn.run(cmd)
        newlist.append('%s/%s' % (transferDir, cmn.lastName(fn)))
    return newlist


def nonGap_char(fn):
    seq = ''.join([line.strip() for line in cmn.file2lines(fn)
        if line.strip() != '' and (not line[0] == '>') ])
    seq = seq.replace('N', '-')
    N = len(seq) - seq.count('-')
    return N


def tell_best_fa(alist):
    #tmp = [each for each in alist
    #        if '/butterfly/' not in each]
    #if len(tmp) != 0:
    #    if len(tmp) == 1:
    #        return tmp[0]
    #    alist = tmp
    print('warnning: telling best fasta from the following sequences')
    print('\n'.join(alist))

    #transfer alea files into here
    alea_files = []
    biohpc_files = []

    for each in alist:
        if '/butterfly/' in each:
            alea_files.append(each)
        else:
            biohpc_files.append(each)

    if alea_files != []:
        new_files = transfer_alea_files(alea_files)
        alist = biohpc_files + new_files

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
        ass = '_'.join(items[1:-4]).replace('_withMito', '').replace('_mitogenome', '').replace('_withNCBImito', '').replace('_NCBImitogenome', '').replace('_with_mito', '')
        if ass == 'cne':
            ass = '3574_assembly_v1'
        try:
            adict[ass] += 1
        except KeyError:
            adict[ass] = 1
    return adict



if __name__=='__main__':
    #options=parse_options()
    try:
        fn=sys.argv[1]
    except:
        print("Usage: *.py NickList requred_keywords", file=sys.stderr)
        sys.exit()

    words = sys.argv[2:]

    IDs = set(cmn.file2lines(fn))

    #note: cne is equal to 3574_assembly_v1
    wdirs = [
            '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/map2fasta',
            '/project/biophysics/Nick_lab/mtang/unbias_SNPs/*/step4_postprocessing/map2fasta'
            ]

    alea_list = cmn.cmd2lines('ssh butterfly@toxea.swmed.edu "ls /archive/butterfly/unbias_pipeline_info/step4_postprocessing/map2fasta/*MITO*.fa"')


    missing = []
    falist = []
    required = ''
    if len(words) != 0:
        required = '|' + '|'.join(['grep %s' % word for word in words])

    faDict = {}
    for ID in IDs:
        taken = []
        for wdir in wdirs:
            cmd = 'ls %s/%s*_m2s.fa  2> /dev/null| grep MITO %s' % (wdir, ID, required)
            taken += cmn.cmd2lines(cmd)

        taken += [each for each in alea_list
                if cmn.lastName(each).split('_')[0] == ID]

        if len(words) != 0:
            taken = [each for each in taken
                    if all([word in each for word in words])]
        faDict[ID] = taken

    #print taken
    all_fa = sum(list(faDict.values()), [])
    ass_count = count_ass_appearance(all_fa)
    best_ass = max(list(ass_count.keys()), key=lambda x: ass_count[x])
    print('the most common assembly is %s, only take fa mapped to this assembly' % best_ass)
    cmn.write_file(best_ass, 'best_assembly.txt')
    sys.exit()

    for ID in faDict:
        alist = faDict[ID]
        taken = [each for each in alist
                if best_ass in each.replace('_withMito', '')]

        #print ID, taken

        if best_ass == 'cne' and len(taken) == 0:
            taken += [each for each in alist
                    if '3574_assembly_v1' in each]

        if best_ass == '3574_assembly_v1' and len(taken) == 0:
            taken += [each for each in alist
                    if 'cne' in each]


        if len(taken) == 0:
            missing.append(ID)
        elif len(taken) == 1:
            falist.append(taken[0])
        else:
            falist.append(tell_best_fa(taken))

    alea_files = []
    biohpc_files = []
    for each in falist:
        if '/archive/butterfly/' in each or ('jshen/h' in each):
            alea_files.append(each)
        else:
            biohpc_files.append(each)

    new_files = transfer_alea_files(alea_files)
    falist = biohpc_files + new_files

    #if len(missing) != 0:
        #try to look for refgenomes
    #    fns = cmn.cmd2lines('ls /work/biophysics/mtang/SNP_calling/indexed_references/mitogenomes/*.fa')
    #    addback = [fn for fn in fns if cmn.lastName(fn).split('_')[0] in missing]
    #    missing = set(missing) - set([cmn.lastName(fn).split('_')[0] for fn in addback])
    #    falist += addback

    if len(missing) != 0:
        print('ATTENTION! the following ID missing sequence!')
        print('\n'.join(missing))
        cmn.write_lines(missing, 'missingMITOs')

    falist.append('')
    cmn.write_lines(falist, 'falist.mito')


    #if len(alea_files) != 0:
    #    print 'the following files need to transfer from /archive server'
    #    print '\n'.join(alea_files)

