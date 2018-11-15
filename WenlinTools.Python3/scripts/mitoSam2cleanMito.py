#!/usr/bin/env python

#function:
#input:
#output:
#algorithm:
#author:wenlin; Date:2012-

import sys
python_lib = '/home/wenlin/my_programs/python_lib'
if python_lib not in sys.path:
    sys.path.append(python_lib)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import cmn
import os
import pysam

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#main
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def find_reference(fn):
    pdir = '/'.join(fn.split('/')[:-1])
    fass_label = '%s/assembly_selfref_v2' % pdir
    fass = '%s/assembly_selfref_v2.fa' % pdir
    if not cmn.filexist(fass):
        print('WARNING: can not find assembly_selfref_v2.fa, use the orginal one')
        reflabel = '_'.join(fn.split('/')[-2].split('_')[1:])
        fass_label = '/work/biophysics/mtang/SNP_calling/indexed_references/%s' % reflabel
    return fass_label



def old_find_reference(fn):
    cmd = 'samtools view -H %s| grep "@PG"| grep bwa' % fn
    info = cmn.cmd2info(cmd)
    items = info.strip().split()
    for i, item in enumerate(items):
        if item == '-M':
            ref = items[i+1]
            break

    if ref[-3:] == '.fa':
        ref = ref[:-3]

    print('found ref: %s' % ref)
    return ref

def find_length_file(splabel):
    items = splabel.split('_')
    reflabel = '_'.join(items[1:])
    flength = '/work/biophysics/mtang/SNP_calling/indexed_references/%s_scafLength.txt' % reflabel
    return flength

def bam2statistics(fn, preLabel, splabel):
    samfile = pysam.AlignmentFile(fn, 'rb')
    total = 0 #total base pair in the bam
    total_count = set([])#total reads
    aligned_bp = [0, 0] #aligned all, 50%reads
    aligned_count = [set([]), set([])] #aligned reads, aligned 50% reads
    for record in samfile:
        ID = record.query_name
        #only take mito
        scaf = record.reference_name
        if 'mito' not in scaf:
            continue

        if record.is_read1:
            ID += '/1'
        elif record.is_read2:
            ID += '/2'

        if record.is_secondary:
            continue

        read_length = record.query_length
        if ID not in total_count:
            total += read_length

        total_count.add(ID)

        if record.is_unmapped:
            continue

        #reach here if they are mapped
        aligned_count[0].add(ID)
        aligned_pairs = record.get_aligned_pairs()
        aligned_length = 0
        for i, j in aligned_pairs:
            if i != None and (j != None):
                aligned_length += 1

        aligned_bp[0] += aligned_length

        if float(aligned_length) / read_length < 0.5:
            continue

        #reach here if more than 50% mapped
        aligned_count[1].add(ID)
        aligned_bp[1] += aligned_length

    #header = '\t'.join('sample label totalN totalBp alignedN alignedBp aligned50%'.split())
    stat = []
    line = [splabel, preLabel, len(total_count), total, len(aligned_count[0]), aligned_bp[0], len(aligned_count[1]), aligned_bp[1]]
    stat.append('\t'.join(map(str, line)))
    return stat


if __name__=='__main__':
    #options=parse_options()
    try:
        fn=os.path.abspath(sys.argv[1])
        fref = os.path.abspath(sys.argv[2])
        outdir = sys.argv[3]
    except:
        print("Usage: *.py realigned.sam fref outlabel", file=sys.stderr)
        sys.exit()

    #/archive/butterfly/SNP_results/debiased/14063C05_1504_mito_rotated/realigned_reads_step2.bam
    #splabel = fn.split('/')[-2]
    splabel = outdir
    wdir = outdir
    cmn.mkdir(wdir)
    os.chdir(wdir)
    cwd = os.getcwd()

    cmn.run('ln -s %s ' % fn)

    fsam = cmn.lastName(fn)
    fin = fsam[:-4] + '.bam'
    cmd = 'module add samtools; samtools view -h -b %s -o tmp.bam; samtools sort tmp.bam -o %s; samtools index %s' % (fsam, fin, fin)
    cmn.run(cmd)

    header = '\t'.join('sample label totalN totalBp alignedN alignedBp aligned50%N aligned50%Bp'.split())
    stat_info = [header]
    stat_info += bam2statistics(fin, 'before filtering', splabel)

    #cmd = 'samtools index %s; samtools view -h %s > %s.sam' % (fin, fin, splabel)
    #cmn.run(cmd)

    #assuming the bam file has been indexed
    #assuming the mito scaffold has a name with '_mito'
    print('preprocessing data for gatk run...')
    fsam = 'mitoOnly.sam'
    fsamHead = 'samHeader.txt'
    cmd = 'samtools view -H %s > %s;cp %s %s; samtools view %s | grep _mito >> %s' % (fin, fsam, fsam, fsamHead, fin, fsam)
    cmn.run(cmd)

    #key step, doing the filtering
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/filter_highQual_mapping.py %s > filter.log' % (fsam)
    cmn.run(cmd)

    cmd = 'mv %s nextStep.sam; cat filtered.sam >> nextStep.sam; rm filtered.sam;' % fsamHead
    cmn.run(cmd)

    #@PG     ID:bwa  VN:0.7.12-r1039 CL:/home2/wli/local/bwa-0.7.12/bwa mem -t 32 -M /project/biophysics/Nick_lab/wli/sequencing/verify_barcodes/Astraptes_barcodeCaseStudy/clean_ref_assembly/3614_assembly_v2sp_withCompMito 15103D03_R1
    cmd = 'samtools view -h -b nextStep.sam -o nextStep.bam; samtools index nextStep.bam'
    cmn.run(cmd)

    #assuming the ref in the same directory as bam and named self something
    #fref_label = find_reference(fn)
    cmn.run('ln -s %s' % fref)
    ref_label = cmn.lastName(fref)[:-3]

    template = cmn.txt_read('/project/biophysics/Nick_lab/wli/sequencing/scripts/templates/re-SNPcall.template')
    info = template.replace('[ref_label]', ref_label)
    info = info.replace('[input.bam]', 'nextStep.bam')
    fvcf = splabel + '.vcf'
    info = info.replace('[out.vcf]', fvcf)
    info = info.replace('[CWD]', cwd)

    #os.chdir('..')
    #cmn.mkdir('job_files')
    #fjob = 'job_files/gc%s.job' % (splabel)
    fjob = 'gc%s.job' % (splabel)
    cmn.write_file(info, fjob)
    #print 'please submit %s to queue' % fjob

    print('running gatk job...')
    cmd = 'bash %s &> SNP.log' % (fjob)
    cmn.run(cmd)

    print('postprocessing gatk data...')
    #then process the vcf file
    flength = find_length_file(splabel)
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/vcf2map_withRefLength_mitoOnly.py %s %s' % (fvcf, flength)
    cmn.run(cmd)

    fmap = '%s_MITOonly.map' % splabel
    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/map2fasta_mito.py %s ' % fmap
    cmn.run(cmd)

    cmd = 'python /work/biophysics/mtang/SNP_calling/scripts/vcf2coverage_mitoOnly.py %s > %s.vcfCov' % (fvcf, splabel)
    cmn.run(cmd)

    stat_info += bam2statistics('realigned_reads_clean.bam', 'after filtering', splabel)
    dn = '%s.samStat' % (splabel)
    cmn.write_lines(stat_info, dn)
    #cmn.write_file('nextStep.sam', 'filelist')
    #cmd = 'python /export/home8/wenlin/fixEudamineMito/confident_calling/step4_parse_sams.py filelist > parse.log'
    #cmn.run(cmd)

    #cmd = 'touch fbait; touch bad_aln; python /export/home8/wenlin/fixEudamineMito/confident_calling/step7_threading_good_reads.py 2nd_sam_aln.txt fbait bad_aln %s > threading.log' % splabel
    #cmn.run(cmd)


