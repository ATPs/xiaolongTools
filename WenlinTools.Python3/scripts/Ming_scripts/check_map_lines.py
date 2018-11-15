import sys
import os


python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

def read_refN(ref_genomes):
    adict = {}
    for ref in ref_genomes:
        fN = '%s/%s_scaf_header.lines' % (ref_dir, ref)
        #print fN
        if not cmn.filexist(fN):
            fhead = '%s/%s_scaf.header' % (ref_dir, ref)
            cmd = 'wc -l %s > %s' % (fhead, fN)
            cmn.run(cmd)
        N = int(cmn.txt_read(fN).split()[0])
        adict[ref] = N    
    return adict        


def detect_ref_genomes(sps, wdirs):
    adict = {}
    rset = set([])
    isbad = False
    for sp in sps:
        tmpFns = ['%s/%s/best_mapping.txt' % (wdir, sp)
                    for wdir in wdirs]
        
        bestFns = [fn for fn in tmpFns
                if os.path.exists(fn)]
        
        if len(bestFns) == 0:
            print('Error! no best reference info found for %s' % sp)
            isbad = True

        elif len(bestFns) != 1:
            print('Error! can not decide the reference for %s' % sp)
            print('Please remove the duplications in ')
            print('\n'.join(bestFns))
            isbad = True

        ref = cmn.txt_read(bestFns[0])
        adict[sp] = ref
        rset.add(ref)
    
    if isbad:
        sys.exit()
    return rset, adict        


#1. read in data
fns = cmn.getid(sys.argv[1])

#bwa_dirs = cmn.getid(sys.argv[2])
bwa_dirs = [line.strip().rstrip('/') for line in cmn.getid(sys.argv[2])]

#2. check which reference they used
#sp is unique
mapF_dict = {}
ref_genomes, refmapping = set([]), {}
for fn in fns:
    fnlabel = cmn.lastName(fn).replace('_snp_step2.map', '').replace('_snp_step2_MITO.map', '')
    sp = fnlabel
    items = fnlabel.split('_')
    ref = '_'.join(items[1:])
    ref_genomes.add(ref)
    refmapping[sp] = ref
    try:
        mapF_dict[sp].append(fn)
    except KeyError:
        mapF_dict[sp] = [fn]

#sps = mapF_dict.keys()

#ref_genomes, refmapping = detect_ref_genomes(sps, bwa_dirs)

#3. make the length check
ref_dir = '/work/biophysics/mtang/SNP_calling/indexed_references'

print('validating map files...')
refNdict = read_refN(ref_genomes)

good_maps = []
bad_maps = []
for sp in refmapping:
    ref = refmapping[sp]
    fmaps = mapF_dict[sp]
    refN = refNdict[ref]
    mapN = 0
    for fmap in fmaps:
        N = int(cmn.cmd2info('wc -l %s' % fmap).split()[0])
        mapN += N
    if refN != mapN:
        print('Error! the line of map doesn\'t agree with reference for %s' % sp)
        print('Nref vs Nmap: %s %s; ref is %s\n' % (refN, mapN, ref))
        bad_maps += fmaps
    else:
        good_maps += fmaps


cmn.write_lines(good_maps, 'good_maps.txt')
cmn.write_lines(bad_maps, 'bad_maps.txt')
