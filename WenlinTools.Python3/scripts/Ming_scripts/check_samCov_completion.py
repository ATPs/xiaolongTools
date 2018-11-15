import sys, os

sys.path.append('/work/biophysics/mtang/SNP_calling/scripts')

import cmn

#wdir = sys.argv[1] # step2_bwa_mapping

#wdir = wdir.rstrip('/')
wdir = '.'

fns = cmn.cmd2lines('ls %s/*/*/*.sam' % wdir)

dirs = set(['/'.join(fn.split('/')[:-2]) for fn in fns])
cov_files = cmn.cmd2lines('ls mapped_reads_count/*_cov.count 2> /dev/null')

finished_dirs = set([cmn.lastName(fn).replace('_cov.count', '')
                    for fn in cov_files])

cwd = os.getcwd()


isGood = True

cmds = ['cd %s' % wdir]
for dir in dirs:
    sp = cmn.lastName(dir)
    if sp in finished_dirs and cmn.filexist('mapped_reads_count/%s_cov.count' % sp):
        continue
    isGood = False    
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/tell_best_mapping.py %s &' % dir
    cmds.append(cmd)


cmds.append('\nwait\n')

outdir = '%s/mapped_reads_count' % wdir
cmn.mkdir(outdir)

for dir in dirs:
    sp = cmn.lastName(dir.rstrip('/'))
    dn = '%s/%s_cov.count' % (outdir, sp)
    if sp in finished_dirs and cmn.filexist(dn):
        continue
    isGood = False    

    #cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/count_sam_aligns_bestDir.py %s > %s &' % (dir, dn)
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/collect_sam_aligns_byDir.py %s > %s &' % (dir, dn)
    cmds.append(cmd)

cmds.append('\nwait\n')

if isGood:
    print('good news! everything looks good!')
else:
    fcmd = 'covStatAdd.cmds'
    cmn.write_lines(cmds, fcmd)

    fjob = 'covStatAdd.job'
    cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -p 256GB -t 20 > %s' % (fcmd, fjob)
    cmn.run(cmd)

    print('Error!!!!!')
    print('There are still missing sam coverage report check!')
    print('please submit the job %s to recheck those missing ones' % fjob)
    print('If you still see this Error after your second trial, Please contact Wenlin immediately!!!')

    
