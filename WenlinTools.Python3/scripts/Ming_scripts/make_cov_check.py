import sys, os

sys.path.append('/work/biophysics/mtang/SNP_calling/scripts')

import cmn

wdir = sys.argv[1] # step2_bwa_mapping

wdir = wdir.rstrip('/')

fns = cmn.cmd2lines('ls %s/*/*/*.sam' % wdir)

dirs = set(['/'.join(fn.split('/')[:-2]) for fn in fns])

cwd = os.getcwd()

outdir = '%s/mapped_reads_count' % wdir
cmn.mkdir(outdir)

cmds = []
for dir in dirs:
    cmd = 'cd %s; ' % wdir
    cmd += 'python /work/biophysics/mtang/SNP_calling/scripts/tell_best_mapping.py %s ; ' % dir

    sp = cmn.lastName(dir.rstrip('/'))
    dn = '%s/%s_cov.count' % (outdir, sp)

    #cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/count_sam_aligns_bestDir.py %s > %s &' % (dir, dn)
    cmd += 'python /work/biophysics/mtang/SNP_calling/scripts/collect_sam_aligns_byDir.py %s > %s ; ' % (dir, dn)
    cmds.append(cmd)

fcmd = 'covStat.cmds'
cmn.write_lines(cmds, fcmd)

#fjob = 'covStat.job'
#cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s -p 256GB -t 20 > %s' % (fcmd, fjob)
#cmn.run(cmd)

print('please submit the cmds %s using submit_jobs.py, the result would be in %s' % (fcmd, outdir))

    
