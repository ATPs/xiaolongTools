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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



if __name__=='__main__':
    #options=parse_options()
    try:
        wdir=sys.argv[1]
    except:
        print("Usage: *.py 2_gatk", file=sys.stderr)
        sys.exit()

    #keep:
    #1. py files
    #2. vcf files
    #3. realigned.bam

    rmfiles = [
            'gatk*.job', '*.err', '*.out', '*/*.log', 'ph*.job',
            '*/*_snp.vcf.idx',
            '*/dedup_reads.bai',
            '*/dedup_reads.bam',
            '*/dedup_reads_fix.bai',
            '*/metrics.txt',
            '*/realigned.pileup',
            '*/sorted_reads.bam',
            '*/target_intervals.list',
            '*/tmp -r',
            '*_files -r',
            '*/assembly_selfref.fa',
            '*.sam',
            '*/*.sam',
            ]

    for label in rmfiles:
        cmd = 'rm %s/%s' % (wdir, label)
        print(cmd)

    #remove the indexes
    cmd = 'ls %s/assembly_*.*| grep -v fa|xargs rm' % wdir
    print(cmd)

    cmd = 'ls %s/*/assembly_selfref_v2.*| grep -v fa|xargs rm' % wdir
    print(cmd)
