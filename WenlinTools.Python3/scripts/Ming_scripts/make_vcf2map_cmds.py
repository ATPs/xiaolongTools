import sys
import os


python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn

def detect_ref_genomes(sps, wdirs):
    adict = {}
    rset = set([])
    isbad = False
    for sp in sps:
        tmpFns = ['%s/%s/best_mapping.txt' % (wdir, sp)
                    for wdir in wdirs]
        bestFns = [fn for fn in tmpFns
                if os.path.exists(fn)]
        
        #print sp, tmpFns, bestFns            
        #print sp
        #print '\t'.join(tmpFns)

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

bwa_dirs = [line.strip().rstrip('/') for line in cmn.getid(sys.argv[2])]

#2. check which reference they used
#sp is unique
vcf_dict = {cmn.lastName(fn).replace('_snp_step2.vcf', ''): fn for fn in fns}
#6188_3842_assembly_v2_snp_step2.vcf
#vcf_dict = {cmn.lastName(fn).replace('_snp_step2.vcf', ''): fn for fn in fns}
sps = list(vcf_dict.keys())

#ref_genomes, refmapping = detect_ref_genomes(sps, bwa_dirs)
ref_genomes, refmapping = set([]), {}
for fn in fns:
    #../../step3_gatk/5729_3614_assembly_v1/5729_3614_assembly_v1_snp_step2.vcf
    fnlabel = cmn.lastName(fn).replace('_snp_step2.vcf', '')
    items = fnlabel.split('_')
    sp = items[0]
    ref = '_'.join(items[1:])
    ref_genomes.add(ref)
    refmapping[fnlabel] = ref

cmn.pickle_write(refmapping, 'ref_mapping.dict.pkl')
info = ['%s\t%s\n' % (sp, refmapping[sp]) for sp in refmapping]
cmn.write_file(''.join(info), 'ref_mapping.txt')

#3. make the length check
ref_dir = '/work/biophysics/mtang/SNP_calling/indexed_references'

unChecks = [ref for ref in ref_genomes
    if not os.path.exists('%s/%s_scafLength.txt' % (ref_dir, ref))
    or (not os.path.exists('%s/%s_scaf.header' % (ref_dir, ref)))]

#unChecks = ref_genomes    

print('#######################################################')

ischeck = True

if len(unChecks) == 0:
    print('good news! all reference genome has been preprocessed')
else:
    
    ischeck = False
    cmds = []
    for ref in unChecks:
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/assembly_scaf_length.py %s/%s.fa &' % (ref_dir, ref)
        cmds.append(cmd)
        
        cmd = '/work/biophysics/mtang/SNP_calling/scripts/build_scaf_header.py %s/%s.fa &' % (ref_dir, ref)
        cmds.append(cmd)
    
    
    cmds.append('\nwait\n')
    fcmds = 'checkL.cmds'
    cmn.write_lines(cmds, fcmds)
    
    fjob = 'checkL.job'
    cmd = '/work/biophysics/mtang/SNP_calling/scripts/decorate_job.py %s > %s' % (fcmds, fjob)
    cmn.run(cmd)
    print('%s genomes needs to check length, please submit %s to queue' % (len(unChecks), fjob))

print('#######################################################')

if not ischeck:
    print('\nPlease re-run this script after finish running checkL.job\n')
    sys.exit()


#4. make the command to get the map files
#also check the line of vcf files to see if the line meets. otherwise something was wrong

#print 'validating vcf files...'
#refNdict = {}
#for ref in ref_genomes:
#    fhead = '%s/%s_scaf.header' % (ref_dir, ref)
#    N = int(cmn.cmd2info('wc -l %s' % fhead).split()[0])
#    refNdict[ref] = N
#    print ref, fhead, N


cmds = []
for sp in refmapping:
    ref = refmapping[sp]
    fref = '%s/%s_scafLength.txt' % (ref_dir, ref)
    fvcf = vcf_dict[sp]
    #refN = refNdict[ref]
    #vcfN = int(cmn.cmd2info('grep -v "^#" %s | wc -l' % fvcf).split()[0])
    #if refN != vcfN:
    #    print 'Error! the line of vcf doesn\'t agree with reference:\n%s ' % fvcf 
    #    print 'Nref vs Nvcf: %s %s; ref is %s' % (refN, vcfN, fref)
    cmd = '/work/biophysics/mtang/SNP_calling/scripts/vcf2map_withRefLength_checkMito.py %s %s' % (fvcf, fref)
    cmds.append(cmd)


cmds.append('')
dn = 'vcf2map.cmds'
cmn.write_lines(cmds, dn)

print('After finish checking scaffold length, please use following command to submit jobs')
print('\n>>> /work/biophysics/mtang/SNP_calling/scripts/submit_jobs.py %s [#node] vcf2map -p 256GB\n' % dn)
print('-p specifies the partition it submitted to')
print('[#node] is the number of nodes and should be adjusted according to number of lines in %s' % dn)

