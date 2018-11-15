import sys
python_lib = '/work/biophysics/mtang/SNP_calling/scripts'
if python_lib not in sys.path:
    sys.path.append(python_lib)

import cmn
import os

def check_fastq_bp(fn):
    count = 0
    with open(fn) as fp:
        for i, line in enumerate(fp):
            if i % 4 == 1:
                count += len(line.strip())
    return count                

def format_readSize(count):
    Gunit = 1000000000
    Munit = 1000000
    Kunit = 1000
    if count > Gunit/100:
        s = '%.2fGbp' % (float(count)/Gunit)
    elif count > Munit/100:
        s = '%.2fMbp' % (float(count)/Munit)
    else:
        s = '%.2fMbp' % (float(count)/Kunit)
    return s

def get_ref_length(reflabel):
    flength = '/work/biophysics/mtang/SNP_calling/indexed_references/%s_scafLength.txt' % reflabel
    total = 0
    for line in cmn.file2lines(flength):
        total += int(line.strip().split()[-1])
    return total

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fq_list = sys.argv[1].split(',')
vcf_fn = sys.argv[2]
samdir = sys.argv[3]
vcfCovDir = sys.argv[4]

refresh = any([each=='-r' for each in sys.argv])
print(refresh)

#1. step1, get the fastq total data size
readSize = 0
if fq_list[0] == 'fake.fq':
    readSizeStr = 'NA'
else:
    for fn in fq_list:
        dataSize = check_fastq_bp(fn)
        readSize += dataSize

    readSizeStr = format_readSize(readSize)
print(readSizeStr)

#2. step2, get the reference and its length
vcf_label = cmn.lastName(vcf_fn).replace('_snp_step2.vcf', '')
items = vcf_label.split('_')
sp = items[0]
reflabel = '_'.join(items[1:])
ref_length = get_ref_length(reflabel)
print(ref_length)

#3. get percentage of mapping
#../../step2_bwa_mapping
#TODO: if sam data available, recompute it
cmd = 'cat %s/mapped_reads_count/*| grep %s| grep %s' % (samdir, sp, reflabel)
info = cmn.cmd2info(cmd)
items = info.strip().split()

if len(items) == 0:
    print('Error! can not find map percentage for %s %s' % (sp, reflabel))
    mapPercentage = 'NA'
else:
    if len(items) == 4:#old format, ignore the mapN
        mapPercentage = 'oldstat'
        mapN, totalN = list(map(int, items[2:4]))
    else:
        mapN, totalN, halfN, pPercent = list(map(float, items[-4:]))
    #mapPercentage = float(mapN) / totalN
    mapPercentage = 'ready'

   
print(mapPercentage)
#print float(readSize) * mapPercentage / ref_length

#4. get vcf stat
dn = '%s/%s_%s_vcf.cov' % (vcfCovDir, sp, reflabel)
if os.path.exists(dn):
    cmd = 'cat %s| grep report' % (dn)
    info = cmn.cmd2info(cmd)
    items = info.strip().split()
else:
    items = []
if len(items) != 6:
    #recompute the cov
    print('coverage info is not sufficent, re-run it and save to %s' % dn)
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/vcf2coverage.py %s > %s' % (vcf_fn, dn)
    cmn.run(cmd)
    items = cmn.file2lines(dn)[-1].strip().split()

vcfPercent, cov, covM = list(map(float, items[3:6]))

header = 'sample reference data_amount map_read%(byReads) map_read%(halfmap) map_read%(byPosition) expected_coverage genome_coverage coverage_mean coverage_median'.split() 
final = ['\t'.join(header)]
line = [sp, reflabel, readSizeStr]
if mapPercentage == 'NA':
    line.append('NA\tNA\tNA\tNA')
else:
    line.append('%.2f%%' % (float(mapN) / totalN * 100) )
    if mapPercentage == 'oldstat':
        line.append('NA\tNA\tNA')
    else:    
        line.append('%.2f%%' % (float(halfN) / totalN * 100) )
        line.append('%.2f%%' % (pPercent * 100) )
        line.append('%.2f' % (float(readSize) * pPercent / ref_length))
line.append('%.2f%%' % (vcfPercent * 100))
line.append('%.2f' % cov)
line.append('%.2f' % covM)
final.append('\t'.join(line))

final.append('')

dn = '%s_%s_stat.report' % (sp, reflabel)
cmn.write_lines(final, dn)




