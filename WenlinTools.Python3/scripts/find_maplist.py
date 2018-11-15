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
        ass = '_'.join(items[1:-3]).replace('_withMito', '')
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
            '/project/biophysics/Nick_lab/mtang/archive/step4_postprocessing/vcf2map',
            '/project/biophysics/Nick_lab/mtang/unbias_SNPs/*/step4_postprocessing/vcf2map'
            ]

    alea_list = cmn.cmd2lines('ssh butterfly@toxea.swmed.edu "ls /archive/butterfly/unbias_pipeline_info/step4_postprocessing/vcf2map/*.map| grep -v mitogenome"')


    missing = []
    falist = []
    required = ''
    if len(words) != 0:
        for word in words:
            word = word.strip()
            if word == '3574_assembly_v1' or word == 'cne':
                word = '"3574_assembly_v1\|cne"'

            required += '| grep %s' % word

    words = set(words)

    faDict = {}
    for ID in IDs:
        taken = []
        for wdir in wdirs:
            cmd = 'ls %s/%s*map  2> /dev/null| grep -v MITO | grep -v mitogenome %s' % (wdir, ID, required)
            taken += cmn.cmd2lines(cmd)

        taken += [each for each in alea_list
                if cmn.lastName(each).split('_')[0] == ID and 'MITO' not in cmn.lastName(each) and 'mitogenome' not in cmn.lastName(each)]

        if len(words) != 0:
            tmp = [each for each in taken
                    if all([word in each for word in words])]
            if ('3574_assembly_v1' in words) and ('cne' not in words):
                tmp += [each for each in taken
                        if 'cne' in each]
            if ('cne' in words) and ('3574_assembly_v1' not in words):
                tmp += [each for each in taken
                        if '3574_assembly_v1' in each]

            taken = tmp

        faDict[ID] = taken

    #print taken
    all_fa = sum(list(faDict.values()), [])
    ass_count = count_ass_appearance(all_fa)
    best_ass = max(list(ass_count.keys()), key=lambda x: ass_count[x])
    print('the most common assembly is %s, only take fa mapped to this assembly' % best_ass)

    for ID in faDict:
        alist = faDict[ID]
        print('alist', alist)
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
        print('checklog', ID, taken)
    alea_files = []
    biohpc_files = []
    for each in falist:
        if '/archive/butterfly/' in each or ('jshen/h' in each):
            alea_files.append(each)
        else:
            biohpc_files.append(each)

    new_files = transfer_alea_files(alea_files)
    falist = biohpc_files + new_files

    falist.append('')
    cmn.write_lines(falist, 'statlist')

    if len(missing) != 0:
        print('ATTENTION! the following ID missing sequence!')
        print('\n'.join(missing))
        cmn.write_lines(missing, 'missingIDs')

    #if len(alea_files) != 0:
    #    print 'the following files need to transfer from /archive server'
    #    print '\n'.join(alea_files)

