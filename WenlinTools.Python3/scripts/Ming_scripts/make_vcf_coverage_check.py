import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

fns = cmn.getid(sys.argv[1])

#1. read in data
spdict = {}
for fn in fns:
    sp = cmn.lastName(fn).replace('_snp_step2.vcf', '')
    try:
        spdict[sp].append(fn)
    except:
        spdict[sp] = [fn]

#2. check for sample duplication

isBad = 0
for sp in spdict:
    fns = spdict[sp]
    if len(fns) > 1:
        print('\nError: found duplicated libs for %s ' % sp)
        print('\n'.join(fns))
        print('\n')
        isBad += 1


if isBad != 0:
    print('\nFound %s duplication, please fix and re-run the code' % isBad)
    sys.exit()
        

#make commands
cmds = []
for sp in spdict:
    fn = spdict[sp][0]
    cmd = '/work/biophysics/mtang/SNP_calling/scripts/vcf2coverage.py %s > %s_vcf.cov' % (fn, sp)
    cmds.append(cmd)

cmds.append('')
dn = 'check.cmds'
cmn.write_lines(cmds, dn)

print('please use following command to submit jobs')
print('\n>>> /work/biophysics/mtang/SNP_calling/scripts/submit_jobs.py %s [#node] vcfCheck -p 256GB\n' % dn)
print('-p specifies the partition it submitted to')
print('[#node] is the number of nodes and should be adjusted according to number of lines in %s' % dn)

