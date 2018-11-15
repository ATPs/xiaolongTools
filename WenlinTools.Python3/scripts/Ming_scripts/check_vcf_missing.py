import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

fns = cmn.getid(sys.argv[1])

wdirs = {}
totalSP = set([])

for fn in fns:
    items = fn.strip().split('/')

    wdir = '/'.join(items[:-2])
    sp = items[-2]
    try:
        wdirs[wdir].append(sp)
    except KeyError:
        wdirs[wdir] = [sp]
    
    if sp in totalSP:
        print('Error! duplicated path found for %s' % sp)
        print('please fix!')
    else:
        totalSP.add(sp)

#check each dir
isGood = True
for wdir in wdirs:
    taken_sp = set(wdirs[wdir])
    missing = []

    cmd = 'ls -d %s/*/ | grep -v job_file' % wdir
    sps = [line.split('/')[-2] for line in cmn.cmd2lines(cmd)]

    for sp in sps:
        if sp not in taken_sp:
            missing.append(sp)
    
    if len(missing) != 0:
        isGood = False
        print('\nError! missing vcf for the following samples:')
        print('scaned diretory: %s' % wdir)
        print('missing samples: %s' % (', '.join(missing)))
        print('\n')

if isGood:
    print('Good news! everything looks good!')
