import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

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
bwa_dirs = cmn.getid(sys.argv[1])
fgood = 'good_maps.txt'
fbad = 'bad_maps.txt'

mito = set([])
genome = set([])

for line in cmn.file2lines(fgood):
    sp = line.split('_')[0]
    if 'MITO' in line:
        mito.add(sp)
    else:
        genome.add(sp)


refs, refdict = detect_ref_genomes( genome , bwa_dirs)

#group by ref
refgroups = {}

for sp in genome:
    ref = refdict[sp]
    try:
        refgroups[ref].append(sp)
    except:
        refgroups[ref] = [sp]

#write out output
new = []
keys = sorted(list(refgroups.keys()), key=lambda x: len(refgroups[x]))
for ref in keys:
    sps = refgroups[ref]
    new.append('The following samples are mapped to %s' % ref)
    for sp in sps:
        new.append(sp)
    new.append('')

new.append('The following samples have mito mapped')
for sp in mito:
    new.append(sp)

new.append('\n')

badsps = set([line.split('_')[0] for line in cmn.file2lines(fbad)])
if len(badsps) != 0:
    new.append('The following IDs are not finished yet')
    for sp in badsps:
        new.append(sp)

print('\n'.join(new))



